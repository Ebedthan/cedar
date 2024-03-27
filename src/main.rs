pub mod cli;
pub mod dist;
pub mod sketch;

use std::env;
use std::fs;
use std::io::{self, Write};

use anyhow::Result;
use finch::{open_sketch_file, serialization::Sketch};
use speedytree::{to_newick, Canonical, NeighborJoiningSolver, RapidBtrees};

fn main() -> Result<()> {
    // Read command-line arguments
    let matches = cli::build_cli().get_matches_from(env::args_os());
    let num_threads: usize = matches
        .get_one::<String>("threads")
        .unwrap()
        .parse()
        .unwrap();

    // Create sketches from sequences
    let sketches_path = sketch::generate_sketch_files(&matches);

    // Read created sketches files
    let mut sketches = Vec::new();
    for path in sketches_path {
        sketches.push(open_sketch_file(path)?);
    }
    let new_sketches = sketches.into_iter().flatten().collect::<Vec<Sketch>>();

    // Compute distance between sketches
    let sketch_distance = dist::compute_distances(new_sketches);

    // Compute matrix
    let matrix = dist::sketches_distance_to_matrix(sketch_distance);

    // Compute tree
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

    if !matches.get_flag("keep") {
        fs::remove_dir_all("darwin_tmp")?;
    } else {
        // Create a PHYLIP file to store distances
        dist::to_phylip(matrix.clone())?;
    }

    Ok(())
}
