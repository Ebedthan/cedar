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
pub mod sketch;

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
    let sketches = sketch::create_and_load_sketches(&args.indir, &args.sketch, kmer_size)?;

    // Build tree
    let tree_algorithm = dist::TreeAlgorithm::from_cli(args.canonical, config.threads);

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
    let sketches = sketch::create_and_load_sketches(&args.indir, &args.sketch, kmer_size)?;

    // Compute pairwise distances, each with an uncertainty estimate
    let estimates =
        dist::compute_distances_with_uncertainty(sketches, args.sketch.size, args.sketch.seed);

    let flagged = estimates
        .iter()
        .filter(|e| e.reliability != dist::Reliability::Reliable)
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
