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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mash::distance::compute_distances;
    use itertools::Itertools;
    use std::io::Read;

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

        let result = to_phylip(&dist, temp_dir_path, true);
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

    fn make_sketch_distance(query: &str, reference: &str, mash_distance: f64) -> SketchDistance {
        SketchDistance {
            containment: 0.0,
            jaccard: 0.0,
            mash_distance,
            common_hashes: 0,
            total_hashes: 0,
            query: query.to_string(),
            reference: reference.to_string(),
        }
    }

    #[test]
    fn test_distance_to_matrix_synthetic() {
        // Doesn't depend on the (currently missing) test/sketches fixtures.
        let distances = vec![
            make_sketch_distance("/data/A.fna", "/data/B.fna", 0.1),
            make_sketch_distance("/data/A.fna", "/data/C.fna", 0.2),
            make_sketch_distance("/data/B.fna", "/data/C.fna", 0.3),
        ];

        let matrix = distance_to_matrix(distances);

        assert_eq!(matrix.names.len(), 3);
        let idx = |name: &str| matrix.names.iter().position(|n| n == name).unwrap();
        let (a, b, c) = (idx("A"), idx("B"), idx("C"));

        // Matrix must be symmetric and match the input distances.
        assert_eq!(matrix.matrix[a][b], 0.1);
        assert_eq!(matrix.matrix[b][a], 0.1);
        assert_eq!(matrix.matrix[a][c], 0.2);
        assert_eq!(matrix.matrix[c][a], 0.2);
        assert_eq!(matrix.matrix[b][c], 0.3);
        assert_eq!(matrix.matrix[c][b], 0.3);

        // Diagonal is always 0 (a genome's distance to itself).
        assert_eq!(matrix.matrix[a][a], 0.0);
        assert_eq!(matrix.matrix[b][b], 0.0);
        assert_eq!(matrix.matrix[c][c], 0.0);
    }

    #[test]
    fn test_distance_to_matrix_empty_input() {
        let matrix = distance_to_matrix(vec![]);
        assert!(matrix.matrix.is_empty());
        assert!(matrix.names.is_empty());
    }

    #[test]
    fn test_to_phylip_append_adds_to_existing_content_rather_than_overwriting() {
        let dist1 = DistanceMatrix {
            matrix: vec![vec![0.0, 0.1], vec![0.1, 0.0]],
            names: vec!["X".to_string(), "Y".to_string()],
        };
        let dist2 = DistanceMatrix {
            matrix: vec![vec![0.0, 0.2], vec![0.2, 0.0]],
            names: vec!["Z".to_string(), "W".to_string()],
        };

        let temp_dir = tempfile::tempdir().unwrap();
        let temp_dir_path = temp_dir.path();

        to_phylip(&dist1, temp_dir_path, true).unwrap();
        to_phylip(&dist2, temp_dir_path, true).unwrap();

        let mut file =
            std::fs::File::open(format!("{}/distance.phylip", temp_dir_path.display())).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();

        // Both writes' headers and rows must be present, in order -- proof
        // append=true actually appends rather than silently overwriting.
        let expected = "2\nX 0 0.1\nY 0.1 0\n2\nZ 0 0.2\nW 0.2 0\n";
        assert_eq!(contents, expected);

        temp_dir.close().unwrap();
    }

    #[test]
    fn test_tree_algorithm_from_cli() {
        match TreeAlgorithm::from_cli(true, 4) {
            TreeAlgorithm::Canonical => {}
            other => panic!("expected Canonical, got {:?}", other),
        }
        match TreeAlgorithm::from_cli(false, 4) {
            TreeAlgorithm::RapidBTrees { num_threads } => assert_eq!(num_threads, 4),
            other => panic!("expected RapidBTrees, got {:?}", other),
        }
    }
}
