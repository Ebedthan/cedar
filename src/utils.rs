use crate::dist;
use std::fs;
use std::io::{self, Write};

pub fn compute_newick_tree(
    matrix: &speedytree::DistanceMatrix,
    is_canonical: bool,
    num_threads: usize,
) -> anyhow::Result<String> {
    if is_canonical {
        let tree =
            speedytree::NeighborJoiningSolver::<speedytree::Canonical>::default(matrix.clone())
                .solve()
                .unwrap();
        Ok(speedytree::to_newick(&tree))
    } else {
        let tree =
            speedytree::NeighborJoiningSolver::<speedytree::RapidBtrees>::default(matrix.clone())
                .set_chunk_size(std::cmp::max(matrix.size() / num_threads, 1))
                .solve()
                .unwrap();
        Ok(speedytree::to_newick(&tree))
    }
}

pub fn output_tree(matches: &clap::ArgMatches, newick: String) -> anyhow::Result<()> {
    if matches.contains_id("output") {
        let path = matches.get_one::<String>("output").unwrap();
        let mut file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)?;
        file.write_all(newick.as_bytes())?;
    } else {
        writeln!(io::stdout(), "{}", newick)?;
    }
    Ok(())
}

pub fn manage_tempdir(
    matches: &clap::ArgMatches,
    matrix: &speedytree::DistanceMatrix,
    tempdir: &str,
) -> anyhow::Result<()> {
    if matches.get_flag("keep") {
        dist::to_phylip(matrix.clone(), tempdir)?;
    } else {
        fs::remove_dir_all(tempdir)?;
    }
    Ok(())
}
