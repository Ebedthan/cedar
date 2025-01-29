// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod cli;
pub mod dist;
pub mod sketch;
pub mod utils;

use rayon::prelude::*;
use std::{env, fs, process};

use anyhow::Context;
use finch::serialization::Sketch;

fn main() -> anyhow::Result<()> {
    // Read command-line arguments
    let matches = cli::build_cli().get_matches_from(env::args_os());
    let filenames: Vec<&str> = matches
        .get_many::<String>("INPUT")
        .unwrap_or_default()
        .map(String::as_str)
        .collect();

    // Ensure at least three files are provided
    if filenames.len() < 3 {
        eprintln!("Error: at least three files should be specified");
        process::exit(1);
    }

    // Parse CLI arguments
    let kmer_size = parse_arg::<u8>(&matches, "kmer")?;
    let sketch_size = parse_arg::<usize>(&matches, "size")?;
    let oversketch = parse_arg::<usize>(&matches, "oversketch")?;
    let seed = parse_arg::<u64>(&matches, "seed")?;
    let num_threads = parse_arg::<usize>(&matches, "threads")?;

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Create temporary directory
    let tempdir = "darwin_tmp";
    fs::create_dir_all(tempdir).context(format!("Could not create temp directory: {}", tempdir))?;

    // Step 1: Create sketches from sequences
    // 1.1. Create sketches using CLI args and sketch::create_sketches function
    let sketches_path = sketch::create_sketches(
        &filenames,
        kmer_size,
        sketch_size,
        oversketch,
        seed,
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
    // 3.1. Compute tree
    let is_canonical = matches.get_flag("canonical");
    let newick: String = utils::compute_newick_tree(&matrix, is_canonical, num_threads)?;

    // 3.2. Output tree
    utils::output_tree(&matches, newick)?;

    // Manage tempdir and tempfiles
    utils::manage_tempdir(&matches, &matrix, tempdir)?;

    Ok(())
}

/// Parses an argument from the CLI and converts it to the specified type
fn parse_arg<T: std::str::FromStr>(matches: &clap::ArgMatches, key: &str) -> anyhow::Result<T>
where
    T::Err: std::fmt::Display,
{
    matches
        .get_one::<String>(key)
        .context(format!("Missing required argument: {}", key))?
        .parse::<T>()
        .map_err(|e| anyhow::anyhow!("Invalid value for {}: {}", key, e))
}
