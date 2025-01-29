// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::{
    collections::HashMap,
    fs::{self},
    io::Write,
    path::{Path, PathBuf},
};

use finch::{
    distance::distance,
    serialization::{Sketch, SketchDistance},
};
use itertools::Itertools;
use speedytree::DistanceMatrix;

/// Compute distance between sketches
pub fn compute_distances(sketches: Vec<Sketch>) -> Vec<SketchDistance> {
    sketches
        .into_iter()
        .combinations_with_replacement(2)
        .filter_map(|pair| {
            let dist = distance(&pair[0], &pair[1], false).ok()?;
            (dist.mash_distance <= 1.0).then_some(dist)
        })
        .collect()
}

/// Computes a distance matrice from a list of sketches distances
pub fn distance_to_matrix(distances: Vec<SketchDistance>) -> DistanceMatrix {
    // SketchDistance contains more data than needed for this task like
    // containment, jaccard distance, etc.
    // So I initialise a Vec to store only the needed data from SketchDistance
    // needed data: query name, reference name, mash distance
    let mut map: HashMap<(String, String), f64> = HashMap::new();

    for distance in distances {
        let query_basename = Path::new(&distance.query)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let ref_basename = Path::new(&distance.reference)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        map.insert((query_basename, ref_basename), distance.mash_distance);
    }

    // Collect unique names
    let unique_names: Vec<String> = map.keys().map(|(q, _)| q.clone()).unique().collect();
    let n = unique_names.len();

    // Initialize N x N matrix with zeros
    let mut matrix = vec![vec![0.0; n]; n];

    // Fill the matrix using precomputed distances
    // Hold on! Here something is happening. First the SketchDistance struct
    // DOES NOT CONTAINS all pairwise distances but only non repeating
    // distances. However, the [speedytree] NJ trees functions uses
    // a N x N matrix. So to create this matrix, I generate all the combinations
    // using a cartesian product and then fill the matrix.
    for (i, j) in (0..n).cartesian_product(0..n) {
        matrix[i][j] = *map
            .get(&(unique_names[i].clone(), unique_names[j].clone()))
            .or_else(|| map.get(&(unique_names[j].clone(), unique_names[i].clone())))
            .unwrap_or(&0.0);
    }

    DistanceMatrix {
        matrix,
        names: unique_names,
    }
}

/// Write a PHYLIP file from a distance matrice
pub fn to_phylip(dist: DistanceMatrix, output: &str) -> anyhow::Result<()> {
    let mut file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(PathBuf::from(output).join("distance.phylip"))?;

    writeln!(file, "{}", dist.names.len())?;

    for (name, row) in dist.names.iter().zip(&dist.matrix) {
        writeln!(file, "{} {}", name, row.iter().format(" "))?;
    }

    Ok(())
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
        assert_eq!(distances.len(), 6);

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
        assert_eq!(matrix.matrix.len(), 3);
        assert_eq!(matrix.matrix[0].len(), 3);
        assert_eq!(matrix.names.len(), 3);
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
        let temp_dir_path = temp_dir.path().to_str().unwrap().to_string();

        let result = to_phylip(dist.clone(), &temp_dir_path);
        assert!(result.is_ok());

        // Verify that the output file is created
        let mut file = std::fs::File::open(format!("{}/distance.phylip", temp_dir_path)).unwrap();
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
