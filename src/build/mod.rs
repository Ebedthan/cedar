// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use anyhow::Context;
use std::fs;
use std::fs::File;
use subprocess::Redirection;

use crate::cli;
use crate::utils;

mod bootstrap;
pub mod dist;
pub mod sketch;

pub fn build_tree_from_genomes(args: &cli::BuildArgs, threads: usize) -> anyhow::Result<()> {
    let genomes = args.genomes.as_ref().unwrap();

    // Validate inputs early
    utils::validate_inputs(&genomes)?;

    // Compute genome statistics in parallel
    let stats = utils::compute_genome_stats(&genomes)?;

    // Check for genome outliers and exit early if found
    utils::check_genome_outliers(&stats, 0.01)?;

    // Determine k-mer size
    let kmer_size = utils::determine_kmer_size(&args, &stats);

    // Set bootstrap replicates
    let reps = args.bootstrap.unwrap_or(100);
    if args.bootstrap.is_some() {
        println!("Boostrap replicated: {}", reps);
    }

    // Create temp directory
    fs::create_dir_all(&args.tempdir)
        .with_context(|| format!("Could not create temp directory: {}", &args.tempdir))?;

    // Create sketches
    let sketches = sketch::create_and_load_sketches(&args, kmer_size)?;

    let tree_algorithm = dist::TreeAlgorithm::from_cli(args.canonical, threads);

    if let Some(bootstrap_reps) = args.bootstrap {
        bootstrap::build_bootstrap_tree(
            sketches,
            bootstrap_reps,
            &tree_algorithm,
            args.output.clone(),
        )
    } else {
        bootstrap::build_single_tree(sketches, &tree_algorithm, &args)
    }
}

pub fn build_tree_from_orthologous_groups(
    args: &cli::BuildArgs,
    threads: usize,
) -> anyhow::Result<()> {
    let orthologous = args.groups.as_ref().unwrap();

    // Validate input early
    utils::validate_inputs(&orthologous)?;

    Ok(())
}

fn run_mafft(infile: &str, outfile: &str) -> anyhow::Result<bool> {
    let run = subprocess::Exec::cmd("linsi")
        .arg(infile)
        .stdout(Redirection::File(File::create(outfile)?))
        .join()?
        .success();
    Ok(run)
}

/*
fn run_trimal(infile: &str, outfile: &str) -> anyhow::Result<()> {
    let run = subprocess::Exec::cmd("trimal").arg(infile)
}*/
