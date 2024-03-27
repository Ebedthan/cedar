pub mod cli;
pub mod dist;
pub mod sketch;

use std::env;

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
    let is_canonical = matches.get_flag("canonical");
    if is_canonical {
        let tree = NeighborJoiningSolver::<Canonical>::default(matrix)
            .solve()
            .unwrap();
        let newick = to_newick(&tree);

        println!("{newick}");
    } else {
        let tree = NeighborJoiningSolver::<RapidBtrees>::default(matrix)
            .set_chunk_size(num_threads)
            .solve()
            .unwrap();
        let newick = to_newick(&tree);

        println!("{newick}");
    }

    Ok(())
}
