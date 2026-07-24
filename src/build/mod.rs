// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs;
use std::path::PathBuf;

use crate::cli;
use crate::utils;

mod bootstrap;
pub mod connectivity;
pub mod dist;
pub mod matrix;
use finch::serialization::Sketch;

use crate::build::matrix::TreeAlgorithm;
use crate::dist::rescue::{compute_distances_with_uncertainty, rescue_kmer_candidates};
use crate::mash::sketch::{create_and_load_sketches, k_computing};
use crate::mash::uncertainty::compute_base_distances_with_uncertainty;
use crate::mash::uncertainty::Reliability;

/// Configuration for phylogenetic analysis
#[derive(Debug, Clone)]
pub struct PhyloConfig {
    pub threads: usize,
    pub bootstrap_reps: usize,
    pub outlier_threshold: f64,
    pub buffer_size: usize,
}

impl Default for PhyloConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            bootstrap_reps: 1000,
            outlier_threshold: 0.01,
            buffer_size: 8192, // 8KB
        }
    }
}

/// Result of `resolve_connectivity_with_smaller_kmer`'s global k-mer
/// search for `--include-div-pairs`.
struct ConnectivityResolution {
    sketches: Vec<Sketch>,
    kmer_size: u8,
    fully_connected: bool,
    remaining_load_bearing: Vec<(String, String)>,
}

/// Search for a smaller dataset-wide k-mer size that resolves connectivity
/// gaps identified by `analyze_connectivity`. Re-sketching the *entire*
/// dataset at each candidate (not just the load-bearing pairs) because
/// Mash distance depends directyl on k, so mixing k across a single distance
/// matrix would break NJ's implicit assumption that every entry is comparable.
///
/// The floor is `k_computing` applied to the LARGEST genome in the dataset:
/// the largest genome is the most permissive one for chance k-mer collisions
/// so it sets the k below which any apparent fix would just be k-mer-space
/// exhaustion, not real homology.
///
/// Stops at the first (largest) k where every genome connects via
/// non-divergent pairs alone. If the floor iss reached without achieving that
/// returns the floor's sketches/estimates as best effort, along with
/// whichever pairs are still load-bearing at that point.
fn resolve_connectivity_with_smaller_kmer(
    indir: &str,
    sketch_args: &cli::SketchArgs,
    current_kmer: u8,
    genome_sizes: &[(String, usize)],
) -> Result<ConnectivityResolution> {
    let largest_genome_size = genome_sizes
        .iter()
        .map(|(_, size)| *size as u32)
        .max()
        .unwrap_or(0);

    let floor = k_computing(largest_genome_size, 0.01);

    let mut candidates = vec![current_kmer];
    candidates.extend(rescue_kmer_candidates(current_kmer, floor));
    let last_index = candidates.len() - 1;

    for (i, &k) in candidates.iter().enumerate() {
        let sketches = create_and_load_sketches(indir, sketch_args.size, sketch_args.seed, k)?;
        let base_estimates = compute_base_distances_with_uncertainty(sketches.clone());
        let edges = dist::connectivity_edges(&base_estimates, &sketches, k);
        let load_bearing = connectivity::analyze_connectivity(&edges);

        if load_bearing.is_empty() {
            return Ok(ConnectivityResolution {
                sketches,
                kmer_size: k,
                fully_connected: true,
                remaining_load_bearing: Vec::new(),
            });
        }

        if i == last_index {
            return Ok(ConnectivityResolution {
                sketches,
                kmer_size: k,
                fully_connected: false,
                remaining_load_bearing: load_bearing,
            });
        }
    }

    unreachable!("candidate list always contains at leat current_kmer")
}

pub fn build_tree_using_mash_distance(args: &cli::BuildArgs, threads: usize) -> Result<()> {
    let config = PhyloConfig {
        threads,
        bootstrap_reps: args.bootstrap.unwrap_or(1000),
        ..Default::default()
    };

    // Validate and collect input files
    let inputs = validate_and_collect_inputs(&args.indir)?;

    // Compute genome statistics in parallel and genome size outliers
    let stats = utils::compute_genome_stats(&inputs)?;
    utils::check_genome_outliers(&stats, config.outlier_threshold)?;

    // Determine k-mer size
    let kmer_size = utils::determine_kmer_size(&args.sketch, &stats)?;

    // Create sketches
    let mut sketches =
        create_and_load_sketches(&args.indir, args.sketch.size, args.sketch.seed, kmer_size)?;

    // Evaluate reliability at this k for every pair
    // and check whether excluding divergent pairs would still leave
    // every genome connected via resolvable data alone.
    let base_estimates = compute_base_distances_with_uncertainty(sketches.clone());
    let edges = dist::connectivity_edges(&base_estimates, &sketches, kmer_size);
    let load_bearing = connectivity::analyze_connectivity(&edges);

    if !load_bearing.is_empty() {
        if args.include_div_pairs {
            println!(
                "{} divergent pair(s) are load-bearing for connectivity; \
                searching for a smaller global k-mer size to resolve them...",
                load_bearing.len()
            );
            let resolution = resolve_connectivity_with_smaller_kmer(
                &args.indir,
                &args.sketch,
                kmer_size,
                &stats,
            )?;

            if resolution.fully_connected {
                println!(
                    "Resolved: switching the whole dataset to k = {} for full connectivity.",
                    resolution.kmer_size
                );
            } else {
                println!(
                    "Could not fully resolve connectivity even at the smallest safe k={}; \
                                     proceeding with {} pair(s) still load-bearing and divergent.",
                    resolution.kmer_size,
                    resolution.remaining_load_bearing.len()
                );
                for (q, r) in &resolution.remaining_load_bearing {
                    println!("  - {} vs {}", q, r);
                }
            }

            sketches = resolution.sketches;
        } else {
            println!(
                "Excluding {} divergent pair(s) that are load-bearing for connectivty: ",
                load_bearing.len()
            );
            for (q, r) in &load_bearing {
                println!("  -  {} vs {}", q, r);
            }
            println!(
                "These pairs are the only thing connecting some genomes to the rest of the \
                dataset via resolvable data. Excluding them may still produce a valid tree \
                (NJ can use the remaining edges), but rerun with --include-div-pairs to \
                attempt resolving them with a smaller k-mer size instead."
            );
        }
    }

    // Build tree
    let tree_algorithm = TreeAlgorithm::from_cli(args.canonical, config.threads);

    if let Some(bootstrap_reps) = args.bootstrap {
        println!("Bootstrap replicates: {}", bootstrap_reps);
        bootstrap::build_bootstrap_tree(
            sketches,
            bootstrap_reps,
            &tree_algorithm,
            args.output.clone(),
        )
    } else if let Some(jacknife_reps) = args.jacknife {
        println!(
            "Jackknife replicates: {} (subsampling {:.0}% of hashes per replicate)",
            jacknife_reps,
            args.jacknife_prop * 100.0
        );
        bootstrap::build_jackknife_tree(
            sketches,
            jacknife_reps,
            args.jacknife_prop,
            &tree_algorithm,
            args.output.clone(),
        )
    } else {
        bootstrap::build_single_tree(sketches, &tree_algorithm, args)
    }
}

/// `cedar dist`: compute pairwise Mash distances annotated with an
/// uncertainty estimate on each one. No tree is built — this is for the
/// "is this pair a solid species-boundary call" moment, not for producing
/// a phylogeny.
pub fn compute_pairwise_distances(args: &cli::DistArgs) -> Result<()> {
    // Validate and collect input files
    let inputs = validate_and_collect_inputs(&args.indir)?;

    // Compute genome statistics and flag genome-size outliers up front —
    // an outlier genome size undermines the k-mer size choice for every
    // pairwise distance computed from it, exactly as it would for `build`.
    let stats = utils::compute_genome_stats(&inputs)?;
    utils::check_genome_outliers(&stats, PhyloConfig::default().outlier_threshold)?;

    // Determine k-mer size: manual override, or the same default heuristic
    // `build` uses. `--target-precision`-driven selection isn't implemented
    // yet, so it's rejected explicitly here rather than silently ignored.
    let kmer_size = utils::determine_kmer_size(&args.sketch, &stats)?;

    // Create sketches
    let sketches =
        create_and_load_sketches(&args.indir, args.sketch.size, args.sketch.seed, kmer_size)?;

    // Compute pairwise distances, each with an uncertainty estimate
    let estimates =
        compute_distances_with_uncertainty(sketches, args.sketch.size, args.sketch.seed);

    let flagged = estimates
        .iter()
        .filter(|e| e.reliability != Reliability::Reliable)
        .count();
    println!(
        "{} pairwise distances computed ({} flagged as borderline, unreliable, or lacking shared hashes)",
        estimates.len(),
        flagged
    );

    utils::output_distance_table(args.output.clone(), &estimates)
}

// Validate input files and collect them efficiently
fn validate_and_collect_inputs(indir: &str) -> Result<Vec<PathBuf>> {
    let entries: Vec<_> = fs::read_dir(indir)
        .with_context(|| format!("Failed to read directory: {}", indir))?
        .collect::<Result<Vec<_>, _>>()?;

    let (valid_files, validation_results): (Vec<_>, Vec<_>) = entries
        .into_par_iter()
        .map(|entry| {
            let path = entry.path();
            let is_valid = utils::is_fasta_format(&path);
            let is_multi = if is_valid {
                utils::is_multi_fasta(&path)
            } else {
                false
            };
            (path, (is_valid, is_multi))
        })
        .unzip();

    // Check for validation errors
    let invalid_files: Vec<_> = valid_files
        .iter()
        .zip(validation_results.iter())
        .filter_map(|(path, (is_valid, _))| if !is_valid { Some(path) } else { None })
        .collect();

    if !invalid_files.is_empty() {
        anyhow::bail!(
            "Input validation error: Only FASTA files are allowed. Invalid files: {:?}",
            invalid_files
        );
    }

    let multi_seq_files: Vec<_> = valid_files
        .iter()
        .zip(validation_results.iter())
        .filter_map(|(path, (_, is_multi))| if *is_multi { Some(path) } else { None })
        .collect();

    if !multi_seq_files.is_empty() {
        anyhow::bail!(
            "Input validation error: Multi-sequence FASTA files detected: {:?}",
            multi_seq_files
        );
    }

    Ok(valid_files)
}
