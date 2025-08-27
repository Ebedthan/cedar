// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod bootstrap;
pub mod cli;
pub mod dist;
pub mod nwk;
pub mod sketch;
pub mod utils;

use clap::Parser;
use std::fs;

use anyhow::{Context, Result};

use crate::{
    bootstrap::{build_bootstrap_tree, build_single_tree},
    cli::Commands,
    dist::TreeAlgorithm,
    sketch::create_and_load_sketches,
    utils::{check_genome_outliers, compute_genome_stats, determine_kmer_size, init_rayon_pool},
};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    init_rayon_pool(cli.threads);

    match cli.command {
        Commands::Build(args) => build_command(args, cli.threads),
        Commands::Compare(args) => compare_command(args),
    }
}

fn build_command(args: cli::BuildArgs, threads: usize) -> Result<()> {
    // Validate inputs early
    utils::validate_inputs(&args.input)?;

    // Compute genome statistics in parallel
    let stats = compute_genome_stats(&args.input)?;

    // Check for genome outliers and exit early if found
    check_genome_outliers(&stats, 0.01)?;

    // Determine k-mer size
    let kmer_size = determine_kmer_size(&args, &stats);

    // Set bootstrap replicates
    let reps = args.bootstrap.unwrap_or(100);
    if args.bootstrap.is_some() {
        println!("Boostrap replicated: {}", reps);
    }

    // Create temp directory
    fs::create_dir_all(&args.tempdir)
        .with_context(|| format!("Could not create temp directory: {}", &args.tempdir))?;

    // Create sketches
    let sketches = create_and_load_sketches(&args, kmer_size)?;

    let tree_algorithm = TreeAlgorithm::from_cli(args.canonical, threads);

    if let Some(bootstrap_reps) = args.bootstrap {
        build_bootstrap_tree(sketches, bootstrap_reps, &tree_algorithm, args.output)
    } else {
        build_single_tree(sketches, &tree_algorithm, &args)
    }
}

fn compare_command(args: cli::CompareArgs) -> anyhow::Result<()> {
    // TODO: implement functionality
    println!("Compare command not yet implemented: {:?}", args);
    Ok(())
}
