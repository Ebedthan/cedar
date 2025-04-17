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

pub fn output_tree(output: Option<String>, newick: String) -> anyhow::Result<()> {
    match output {
        Some(path) => {
            let mut file = fs::File::create(path)?;
            file.write_all(newick.as_bytes())?;
        }
        None => {
            writeln!(io::stdout(), "{}", newick)?;
        }
    }
    Ok(())
}

pub fn manage_tempdir(
    keep: bool,
    matrix: &speedytree::DistanceMatrix,
    tempdir: &str,
) -> anyhow::Result<()> {
    if keep {
        dist::to_phylip(matrix.clone(), tempdir)
    } else {
        fs::remove_dir_all(tempdir)?;
        Ok(())
    }
}
