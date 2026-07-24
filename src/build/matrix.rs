use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::mash::distance::extract_basename;
use finch::serialization::SketchDistance;
use speedytree::DistanceMatrix;

/// Computes a distance matrice from a list of sketches distances
pub fn distance_to_matrix(distances: Vec<SketchDistance>) -> DistanceMatrix {
    if distances.is_empty() {
        return DistanceMatrix {
            matrix: vec![],
            names: vec![],
        };
    }

    // Step 1: extract and deduplicate names efficiently
    let mut name_to_index: HashMap<String, usize> = HashMap::new();
    let mut names = Vec::new();

    for dist in &distances {
        let query_basename = extract_basename(&dist.query);
        let ref_basename = extract_basename(&dist.reference);

        for basename in [query_basename, ref_basename] {
            if !name_to_index.contains_key(&basename) {
                let index = names.len();
                name_to_index.insert(basename.clone(), index);
                names.push(basename);
            }
        }
    }

    let n = names.len();
    let mut matrix = vec![vec![0.0; n]; n];

    // step 2: fill matrix using index lookups
    // Hope to have found a better way to fill the matrix which should also be much faster
    for dist in distances {
        let query_basename = extract_basename(&dist.query);
        let ref_basename = extract_basename(&dist.reference);

        if let (Some(&i), Some(&j)) = (
            name_to_index.get(&query_basename),
            name_to_index.get(&ref_basename),
        ) {
            // create the symmetric (n x n) matrix
            matrix[i][j] = dist.mash_distance;
            matrix[j][i] = dist.mash_distance;
        }
    }

    DistanceMatrix { matrix, names }
}

/// Write a PHYLIP file from a distance matrice
pub fn to_phylip(dist: &DistanceMatrix, output: &Path, append: bool) -> anyhow::Result<()> {
    let phylip_path = output.join("distance.phylip");

    let file = if append {
        fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&phylip_path)?
    } else {
        File::create(&phylip_path)?
    };

    let mut writer = BufWriter::new(file);

    // write header
    writeln!(writer, "{}", dist.names.len())?;

    // pre-allocate string buffer and reuse
    let mut line_buffer = String::with_capacity(dist.names.len() * 12);

    for (name, row) in dist.names.iter().zip(&dist.matrix) {
        line_buffer.clear();
        line_buffer.push_str(name);
        line_buffer.push(' ');
        for (i, &value) in row.iter().enumerate() {
            if i > 0 {
                line_buffer.push(' ');
            }
            line_buffer.push_str(&format!("{}", value));
        }
        writeln!(writer, "{}", line_buffer)?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Debug, Clone)]
pub enum TreeAlgorithm {
    Canonical,
    RapidBTrees { num_threads: usize },
}

impl TreeAlgorithm {
    pub fn from_cli(is_canonical: bool, num_threads: usize) -> Self {
        if is_canonical {
            TreeAlgorithm::Canonical
        } else {
            TreeAlgorithm::RapidBTrees { num_threads }
        }
    }
}

pub trait ComputeTree {
    fn compute_newick_tree(&self, algorithm: &TreeAlgorithm) -> anyhow::Result<String>;
}

impl ComputeTree for DistanceMatrix {
    fn compute_newick_tree(&self, algorithm: &TreeAlgorithm) -> anyhow::Result<String> {
        match algorithm {
            TreeAlgorithm::Canonical => {
                let tree = speedytree::NeighborJoiningSolver::<speedytree::Canonical>::default(
                    self.clone(),
                )
                .solve()
                .map_err(|e| anyhow::anyhow!("Tree computation failed: {}", e))?;
                Ok(speedytree::to_newick(&tree))
            }

            TreeAlgorithm::RapidBTrees { num_threads } => {
                let chunk_size = std::cmp::max(self.size() / num_threads, 1);
                let tree = speedytree::NeighborJoiningSolver::<speedytree::RapidBtrees>::default(
                    self.clone(),
                )
                .set_chunk_size(chunk_size)
                .solve()
                .map_err(|e| anyhow::anyhow!("Tree computation failed: {}", e))?;
                Ok(speedytree::to_newick(&tree))
            }
        }
    }
}
