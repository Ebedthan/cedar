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
    let filenames: Vec<&str> = cli.input.iter().map(String::as_str).collect();

    // Ensure all provided files are in fasta format
    if filenames.iter().any(|f| !utils::is_fasta_format(f)) {
        eprintln!("Error: only fasta files are allowed");
        process::exit(1);
    }

    // Check if all files are only single-sequence fasta
    if filenames.iter().any(|f| utils::is_multi_fasta(f)) {
        eprintln!("Error: only single-sequence fasta files are allowed");
        process::exit(1);
    }

    // Ensure at least three files are provided
    if filenames.len() < 3 {
        eprintln!("Error: at least three files should be specified");
        process::exit(1);
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
        &filenames,
        cli.kmer,
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
