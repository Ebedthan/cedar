use crate::dist;
use std::fs::{self, File};
use std::io::BufRead;
use std::io::BufReader;
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

pub fn is_fasta_format(path: &str) -> bool {
    let file = File::open(path).expect("file should exists before opening.");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    while let Some(Ok(line)) = lines.next() {
        let trimmed_line = line.trim();
        if trimmed_line.is_empty() {
            continue;
        }
        return match trimmed_line.chars().next() {
            Some('>') => true,
            _ => false,
        };
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_fasta_format_ok() {
        assert!(is_fasta_format("test/bac168.fna"));
    }

    #[test]
    fn test_is_fasta_format_not_ok() {
        assert!(!is_fasta_format("test/test.fq"));
    }
}
