// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use anyhow::Result;
use std::path::PathBuf;

use crate::cli;
use crate::utils;

mod bootstrap;
pub mod connectivity;
pub mod matrix;
use finch::serialization::Sketch;

use crate::build::connectivity::connectivity_edges;
use crate::build::matrix::TreeAlgorithm;
use crate::dist::rescue::rescue_kmer_candidates;
use crate::mash::sketch::{create_and_load_sketches, k_computing};
use crate::mash::uncertainty::compute_base_distances_with_uncertainty;

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
            outlier_threshold: utils::DEFAULT_OUTLIER_THRESHOLD,
            buffer_size: 8192,
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
/// gaps identified by `analyze_connectivity`. Re-sketches the *entire*
/// dataset at each candidate (not just the load-bearing pairs) because
/// Mash distance depends directly on k, so mixing k across a single
/// distance matrix would break NJ's implicit assumption that every entry
/// is comparable.
///
/// The floor is `k_computing` applied to the LARGEST genome in the
/// dataset: the largest genome is the most permissive one for chance
/// k-mer collisions, so it sets the k below which any apparent fix would
/// just be k-mer-space exhaustion, not real homology.
///
/// Stops at the first (largest) k where every genome connects via
/// non-divergent pairs alone. If the floor is reached without achieving
/// that, returns the floor's sketches/estimates as best effort, along
/// with whichever pairs are still load-bearing at that point.
#[allow(clippy::too_many_arguments)]
fn resolve_connectivity_with_smaller_kmer(
    inputs: &[PathBuf],
    sketch_args: &cli::SketchArgs,
    current_kmer: u8,
    genome_sizes: &[(String, usize)],
    confidence: f64,
    rescue_pvalue: f64,
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
        let sketches = create_and_load_sketches(inputs, sketch_args.size, sketch_args.seed, k)?;
        let base_estimates = compute_base_distances_with_uncertainty(sketches.clone(), confidence);
        let edges = connectivity_edges(&base_estimates, &sketches, k, rescue_pvalue);
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

    unreachable!("candidate list always contains at least current_kmer")
}

pub fn build_tree_using_mash_distance(args: &cli::BuildArgs, threads: usize) -> Result<()> {
    if !(args.uncertainty.confidence > 0.5 && args.uncertainty.confidence < 1.0) {
        anyhow::bail!(
            "--confidence must be strictly between 0.5 and 1.0, got {}",
            args.uncertainty.confidence
        );
    }
    if !(args.uncertainty.rescue_pvalue > 0.0 && args.uncertainty.rescue_pvalue < 1.0) {
        anyhow::bail!(
            "--rescue-pvalue must be strictly between 0.0 and 1.0, got {}",
            args.uncertainty.rescue_pvalue
        );
    }

    let config = PhyloConfig {
        threads,
        bootstrap_reps: args.bootstrap.unwrap_or(1000),
        ..Default::default()
    };

    // Validate and collect input files (files and/or directories)
    let inputs = utils::validate_and_collect_inputs(&args.inputs)?;

    // Compute genome statistics in parallel and genome size outliers
    let stats = utils::compute_genome_stats(&inputs)?;

    let _outliers = utils::check_genome_outliers(&stats, utils::DEFAULT_OUTLIER_THRESHOLD);

    // Determine k-mer size
    let kmer = utils::determine_kmer_size(&args.sketch, &stats)?;

    // Create sketches
    let mut sketches =
        create_and_load_sketches(&inputs, args.sketch.size, args.sketch.seed, kmer.kmer_size)?;

    // Evaluate reliability at this k for every pair and check whether
    // excluding divergent pairs would still leave every genome connected
    // via resolvable data alone.
    let base_estimates =
        compute_base_distances_with_uncertainty(sketches.clone(), args.uncertainty.confidence);
    let edges = connectivity_edges(
        &base_estimates,
        &sketches,
        kmer.kmer_size,
        args.uncertainty.rescue_pvalue,
    );
    let load_bearing = connectivity::analyze_connectivity(&edges);

    if !load_bearing.is_empty() {
        if args.include_div_pairs {
            println!(
                "{} divergent pair(s) are load-bearing for connectivity; \
                 searching for a smaller global k-mer size to resolve them...",
                load_bearing.len()
            );
            let resolution = resolve_connectivity_with_smaller_kmer(
                &inputs,
                &args.sketch,
                kmer.kmer_size,
                &stats,
                args.uncertainty.confidence,
                args.uncertainty.rescue_pvalue,
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
                "Warning: {} divergent pair(s) are load-bearing for connectivity and will \
                 still be used in the tree (NJ requires a complete distance matrix, so these \
                 pairs can't simply be dropped). Treat branches resting on them with caution:",
                load_bearing.len()
            );
            for (q, r) in &load_bearing {
                println!("  - {} vs {}", q, r);
            }
            println!(
                "Rerun with --include-div-pairs to attempt resolving them with a smaller \
                 k-mer size instead."
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

#[cfg(test)]
mod tests {
    use super::*;
    use finch::filtering::FilterParams;
    use finch::sketch_schemes::{KmerCount, SketchParams};

    #[test]
    fn test_phylo_config_defaults() {
        let config = PhyloConfig::default();
        assert_eq!(config.threads, 1);
        assert_eq!(config.bootstrap_reps, 1000);
        assert_eq!(config.outlier_threshold, utils::DEFAULT_OUTLIER_THRESHOLD);
        assert_eq!(config.buffer_size, 8192);
    }

    /// Hand-builds a sketch with a controlled overlap structure: `shared_range`
    /// hashes are common ground for building known Jaccard values against
    /// other sketches, plus `unique_count` hashes starting at `unique_offset`
    /// that are exclusive to this genome. Verified against real
    /// `finch::distance::distance()` before writing this: with sketches all
    /// the same size, this construction reliably produces
    /// jaccard = |shared_range| / sketch_size for any two genomes sharing
    /// the same `shared_range`, as long as each genome's `unique_offset` is
    /// far enough apart that unique ranges never collide or interleave.
    fn make_sketch(
        name: &str,
        shared_range: std::ops::Range<u64>,
        unique_offset: u64,
        unique_count: u64,
        seq_length: u64,
    ) -> Sketch {
        let mut hashes: Vec<KmerCount> = shared_range
            .map(|h| KmerCount {
                hash: h,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            })
            .collect();
        hashes.extend((0..unique_count).map(|i| KmerCount {
            hash: unique_offset + i,
            kmer: vec![],
            count: 1,
            extra_count: 0,
            label: None,
        }));
        hashes.sort_by_key(|k| k.hash); // sortedness is a documented invariant of finch's raw_distance

        Sketch {
            name: name.to_string(),
            seq_length,
            num_valid_kmers: hashes.len() as u64,
            comment: String::new(),
            hashes,
            filter_params: FilterParams::default(),
            sketch_params: SketchParams::Mash {
                kmers_to_sketch: 1000,
                final_size: 1000,
                no_strict: false,
                kmer_length: 21,
                hash_seed: 42,
            },
        }
    }

    #[test]
    fn test_connectivity_pipeline_flags_exactly_one_isolated_genome() {
        // Four genomes (A,B,C,D) mutually share 950/1000 hashes (J=0.95,
        // tightly Reliable). A fifth genome E shares only 5/1000 with each
        // of them (J=0.005, Unreliable) -- E is only reachable through
        // divergent pairs. This exercises the real pipeline
        // (compute_base_distances_with_uncertainty -> connectivity_edges
        // -> analyze_connectivity) end to end on real Sketch/finch
        // distance computation, not mocked Edge/DistanceEstimate values.
        let sketches = vec![
            make_sketch("A", 0..950, 1_000_000, 50, 5_000_000),
            make_sketch("B", 0..950, 2_000_000, 50, 5_000_000),
            make_sketch("C", 0..950, 3_000_000, 50, 5_000_000),
            make_sketch("D", 0..950, 4_000_000, 50, 5_000_000),
            make_sketch("E", 0..5, 9_000_000, 995, 5_000_000),
        ];

        let base_estimates = compute_base_distances_with_uncertainty(sketches.clone(), 0.95);
        let edges = connectivity_edges(&base_estimates, &sketches, 21, 0.01);
        let load_bearing = connectivity::analyze_connectivity(&edges);

        assert_eq!(
            load_bearing.len(),
            1,
            "expected exactly one load-bearing bridge for E"
        );
        assert!(
            load_bearing[0].0 == "E" || load_bearing[0].1 == "E",
            "expected E in the load-bearing pair, got {:?}",
            load_bearing[0]
        );
    }

    #[test]
    fn test_connectivity_pipeline_finds_no_load_bearing_pairs_when_fully_resolvable() {
        // Same four mutually-resolvable genomes, no isolated fifth genome:
        // connectivity should need zero divergent pairs.
        let sketches = vec![
            make_sketch("A", 0..950, 1_000_000, 50, 5_000_000),
            make_sketch("B", 0..950, 2_000_000, 50, 5_000_000),
            make_sketch("C", 0..950, 3_000_000, 50, 5_000_000),
            make_sketch("D", 0..950, 4_000_000, 50, 5_000_000),
        ];

        let base_estimates = compute_base_distances_with_uncertainty(sketches.clone(), 0.95);
        let edges = connectivity_edges(&base_estimates, &sketches, 21, 0.01);
        let load_bearing = connectivity::analyze_connectivity(&edges);

        assert!(load_bearing.is_empty());
    }
}
