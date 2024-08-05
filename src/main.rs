// Copyright 2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod cli;
pub mod dist;
pub mod sketch;
pub mod utils;

use rayon::prelude::*;
use std::env;
use std::fs;

use anyhow::Context;
use finch::serialization::Sketch;

fn main() -> anyhow::Result<()> {
    // Read command-line arguments
    let matches = cli::build_cli().get_matches_from(env::args_os());
    let filenames: Vec<&str> = matches
        .get_many::<String>("INPUT")
        .unwrap_or_default()
        .map(|v| v.as_str())
        .collect();
    if filenames.len() < 3 {
        eprintln!("Error: at least three files should be specified");
        std::process::exit(1);
    }
    let kmer_size: u8 = matches.get_one::<String>("kmer").unwrap().parse()?;
    let sketch_size: usize = matches.get_one::<String>("size").unwrap().parse()?;
    let oversketch: usize = matches.get_one::<String>("oversketch").unwrap().parse()?;
    let seed: u64 = matches.get_one::<String>("seed").unwrap().parse()?;
    let num_threads: usize = matches.get_one::<String>("threads").unwrap().parse()?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Manage temporary directory
    let tempdir = "darwin_tmp";
    fs::create_dir_all(tempdir).context(format!("Could not create tmpdir {}", tempdir))?;

    // Step 1: Create sketches from sequences
    // 1.1. Create sketches using CLI args and sketch::create_sketches function
    let sketches_path =
        sketch::create_sketches(filenames, kmer_size, sketch_size, oversketch, seed, tempdir)?;

    // 1.2. Read created sketches files in a list
    let sketches: Vec<Sketch> = sketches_path
        .into_par_iter()
        .flat_map(|path| finch::open_sketch_file(path).unwrap())
        .collect();

    // Step 2: Compute distance between sketches
    let sketch_distance = dist::compute_distances(sketches);

    // 2.1. Compute matrix
    let matrix = dist::distance_to_matrix(sketch_distance);

    // Step3: Compute tree
    // 3.1. Compute tree
    let is_canonical = matches.get_flag("canonical");
    let newick: String = utils::compute_newick_tree(&matrix, is_canonical, num_threads)?;

    // 3.2. Output tree
    utils::output_tree(&matches, newick)?;

    // Manage tempdir and tempfiles
    utils::manage_tempdir(&matches, &matrix, tempdir)?;

    Ok(())
}
