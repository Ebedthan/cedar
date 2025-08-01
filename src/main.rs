// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod cli;
pub mod dist;
pub mod sketch;
pub mod utils;
use clap::Parser;

use rayon::prelude::*;
use std::{fs, process};

use anyhow::Context;
use finch::serialization::Sketch;

fn main() -> anyhow::Result<()> {
    // Read command-line arguments
    let cli = cli::Cli::parse();
    let filenames = cli.input.as_slice();

    // Validate inputs
    if let Err(e) = utils::validate_inputs(&cli.input) {
        eprintln!("Input validation error: {e}");
        process::exit(1);
    }

    let stats: Vec<(String, usize)> = filenames
        .iter()
        .map(|f| utils::get_seq_stats(f))
        .collect::<anyhow::Result<_>>()?;

    for stat in &stats {
        println!(
            "Genome: {}, size: {}",
            stat.0,
            utils::format_genome_size(stat.1)
        );
    }
    let outliers = utils::detect_outliers(&stats, 0.05)?;
    if !outliers.is_empty() {
        eprintln!("Error: outliers detected in genome sizes");
        for outlier in outliers {
            eprintln!(
                "Genome {} with size {} negatively influence k selection with influential size",
                outlier.0,
                utils::format_genome_size(outlier.1)
            );
        }
        process::exit(1);
    }
    let mut kmer_size = 0_u8;
    if let Some(km) = cli.kmer {
        println!("User-defined k-mer size: {}", km);
    } else {
        let mean_genome_size = &stats.iter().map(|x| x.1 as u32).sum::<u32>() / stats.len() as u32;
        kmer_size = sketch::k_computing(mean_genome_size, 0.01);
        println!(
            "Computed k-mer size (with mean genome size: {} and probability: {}): {}",
            utils::format_genome_size(mean_genome_size as usize),
            0.01,
            kmer_size
        );
    }

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap();

    // Create temporary directory
    let tempdir = "darwin_tmp";
    fs::create_dir_all(tempdir).context(format!("Could not create temp directory: {}", tempdir))?;

    // Step 1: Create sketches from sequences
    // 1.1. Create sketches using CLI args and sketch::create_sketches function
    let sketches_path = sketch::create_sketches(
        filenames,
        kmer_size,
        cli.size,
        cli.oversketch,
        cli.seed,
        tempdir,
    )?;

    // 1.2. Read created sketches files in a list
    let sketches: Vec<Sketch> = sketches_path
        .into_par_iter()
        .flat_map(|path| finch::open_sketch_file(path).unwrap())
        .collect();

    // Step 2: Compute distance between sketches
    let sketch_distance = dist::compute_distances(sketches);

    // 2.1. Compute matrix
    let matrix = dist::distance_to_matrix(sketch_distance);

    // Step 3: Compute tree
    // 3.1. Compute tree;
    let newick: String = utils::compute_newick_tree(&matrix, cli.canonical, cli.threads)?;

    // 3.2. Output tree
    utils::output_tree(cli.output, newick)?;

    // Manage tempdir and tempfiles
    utils::manage_tempdir(cli.keep, &matrix, tempdir)?;

    Ok(())
}
