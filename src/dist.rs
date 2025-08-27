// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

use finch::{
    distance::distance,
    serialization::{Sketch, SketchDistance},
};
use itertools::Itertools;
use speedytree::DistanceMatrix;

/// Compute distance between sketches
/// Uses rayon for parallel processing
pub fn compute_distances(sketches: Vec<Sketch>) -> Vec<SketchDistance> {
    use rayon::prelude::*;
    sketches
        .into_iter()
        .combinations_with_replacement(2)
        .collect::<Vec<_>>() // collect combinations first
        .into_par_iter() // then parallelize distance computation
        .filter_map(|pair| {
            let dist = distance(&pair[0], &pair[1], false).ok()?;
            // early validation to avoid storing invalid distance
            (dist.mash_distance <= 1.0 && dist.mash_distance >= 0.0).then_some(dist)
        })
        .collect()
}

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

/// Cached basename extraction
fn extract_basename(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_string()
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
            line_buffer.push_str(&format!("{:.6}", value));
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile;

    // Test compute_distances function
    #[test]
    fn test_compute_distances() {
        let mut sketches = Vec::new();
        for file in fs::read_dir("test/sketches").unwrap() {
            sketches.push(finch::open_sketch_file(file.unwrap().path()).unwrap());
        }
        let distances = compute_distances(sketches.into_iter().flatten().collect_vec());

        // Assert that the number of distances is correct
        assert_eq!(distances.len(), 3);

        // Assert that each distance is computed correctly
        for distance in &distances {
            assert!(distance.mash_distance <= 1.0);
        }
    }
    // Test distance_to_matrix function
    #[test]
    fn test_distance_to_matrix() {
        let mut sketches = Vec::new();
        for file in fs::read_dir("test/sketches").unwrap() {
            sketches.push(finch::open_sketch_file(file.unwrap().path()).unwrap());
        }
        let data = sketches.into_iter().flatten().collect_vec();

        let distances = compute_distances(data);

        let matrix = distance_to_matrix(distances);

        // Assert that the matrix is computed correctly
        assert_eq!(matrix.matrix.len(), 2);
        assert_eq!(matrix.matrix[0].len(), 2);
        assert_eq!(matrix.names.len(), 2);
    }

    // Test to_phylip function
    #[test]
    fn test_to_phylip() {
        let dist = DistanceMatrix {
            matrix: vec![
                vec![0.0, 0.5, 0.8],
                vec![0.5, 0.0, 0.9],
                vec![0.8, 0.9, 0.0],
            ],
            names: vec![
                "Sketch1".to_string(),
                "Sketch2".to_string(),
                "Sketch3".to_string(),
            ],
        };

        // Create a temporary directory for testing
        let temp_dir = tempfile::tempdir().unwrap();
        let temp_dir_path = temp_dir.path();

        let result = to_phylip(&dist, &temp_dir_path, true);
        assert!(result.is_ok());

        // Verify that the output file is created
        let mut file =
            std::fs::File::open(format!("{}/distance.phylip", temp_dir_path.display())).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();

        let expected_content = "\
            3\n\
            Sketch1 0 0.5 0.8\n\
            Sketch2 0.5 0 0.9\n\
            Sketch3 0.8 0.9 0\n\
        ";
        assert_eq!(contents, expected_content);

        // Clean up the temporary directory
        temp_dir.close().unwrap();
    }
}
