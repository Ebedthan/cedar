// Copyright 2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod cli;
pub mod dist;
pub mod sketch;

use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

use anyhow::format_err;
use finch::{open_sketch_file, serialization::Sketch};
use speedytree::{to_newick, Canonical, NeighborJoiningSolver, RapidBtrees};

fn main() -> anyhow::Result<()> {
    // Read command-line arguments
    let matches = cli::build_cli().get_matches_from(env::args_os());
    let filenames = matches
        .get_many::<String>("INPUT")
        .unwrap_or_default()
        .map(|v| v.as_str())
        .collect::<Vec<_>>();
    let kmer_size: u8 = matches.get_one::<String>("kmer").unwrap().parse()?;
    let sketch_size: usize = matches.get_one::<String>("size").unwrap().parse()?;
    let oversketch: usize = matches.get_one::<String>("oversketch").unwrap().parse()?;
    let seed: u64 = matches.get_one::<String>("seed").unwrap().parse()?;
    let num_threads: usize = matches.get_one::<String>("threads").unwrap().parse()?;

    // Manage temporary directory
    let tempdir = "darwin_tmp";
    if !Path::new(tempdir).exists() {
        fs::create_dir(tempdir).map_err(|_| format_err!("Could not create tmpdir {}", tempdir))?;
    }

    // Step 1: Create sketches from sequences
    // 1.1. Create sketches using CLI args and sketch::create_sketches function
    let sketches_path =
        sketch::create_sketches(filenames, kmer_size, sketch_size, oversketch, seed, tempdir)?;

    // 1.2. Read created sketches files in a list
    let mut sketches = Vec::new();
    for path in sketches_path {
        sketches.push(open_sketch_file(path)?);
    }
    let new_sketches = sketches.into_iter().flatten().collect::<Vec<Sketch>>();

    // Step 2: Compute distance between sketches
    let sketch_distance = dist::compute_distances(new_sketches);

    // 2.1. Compute matrix
    let matrix = dist::distance_to_matrix(sketch_distance);

    // Step3: Compute tree
    // 3.1. Compute tree
    let write_to_stdin = matches.contains_id("output");
    let is_canonical = matches.get_flag("canonical");
    let newick: String = if is_canonical {
        let tree = NeighborJoiningSolver::<Canonical>::default(matrix.clone())
            .solve()
            .unwrap();
        to_newick(&tree)
    } else {
        let tree = NeighborJoiningSolver::<RapidBtrees>::default(matrix.clone())
            .set_chunk_size(num_threads)
            .solve()
            .unwrap();
        to_newick(&tree)
    };

    // 3.2. Output tree
    if !write_to_stdin {
        writeln!(io::stdout(), "{newick}")?;
    } else {
        let path = matches.get_one::<String>("output").unwrap();
        let mut file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)?;
        file.write_all(newick.as_bytes())?;
    }

    // Manage tempdir and tempfiles
    if !matches.get_flag("keep") {
        fs::remove_dir_all(tempdir)?;
    } else {
        // Create a PHYLIP file to store distances
        dist::to_phylip(matrix.clone(), tempdir)?;
    }

    Ok(())
}
