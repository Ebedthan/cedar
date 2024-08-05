// Copyright 2024 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::path::Path;
use std::{collections::HashMap, path::PathBuf};

use finch::{
    distance::distance,
    serialization::{Sketch, SketchDistance},
};
use itertools::Itertools;
use speedytree::DistanceMatrix;
use std::fs;
use std::io::Write;

/// Compute distance between sketches
pub fn compute_distances(sketches: Vec<Sketch>) -> Vec<SketchDistance> {
    let mut distances = Vec::new();
    for skecth_combination in sketches.into_iter().combinations_with_replacement(2) {
        let distance = distance(&skecth_combination[0], &skecth_combination[1], false).unwrap();
        if distance.mash_distance <= 1f64 {
            distances.push(distance);
        }
    }
    distances
}

/// Computes a distance matrice from a list of sketches distances
pub fn distance_to_matrix(distances: Vec<SketchDistance>) -> DistanceMatrix {
    // SketchDistance contains more data than needed for this task like
    // containment, jaccard distance, etc.
    // So I initialise a Vec to store only the needed data from SketchDistance
    // needed data: query name, reference name, mash distance
    let mut needed_data: Vec<(String, String, f64)> = Vec::new();
    for distance in distances {
        needed_data.push((distance.query, distance.reference, distance.mash_distance));
    }
    // Now I store the data grabed above into a HashMap to ease the distance query
    // The query and reference names can come as a path so I always extract the
    // basename to have a nice naming in the tree file
    let mut map: HashMap<(String, String), f64> = HashMap::new();
    for (query, reference, dist) in needed_data.iter() {
        let query_basename = Path::new(&query)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        let ref_basename = Path::new(&reference)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        map.insert((query_basename, ref_basename), *dist);
    }
    // Get only the uniques names
    let mut unique_names: Vec<String> = Vec::new();
    for (query, _) in map.keys() {
        if !unique_names.contains(query) {
            unique_names.push(query.clone());
        }
    }
    let n = unique_names.len();
    // Initialise the matrix with zeros so I can do things like mat[i][j]
    let mut matrix = vec![vec![0.0; n]; n];
    // Hold on! Here something is happening. First the SketchDistance struct
    // DOES NOT CONTAINS all pairwise distances but only non repeating
    // distances. However, the [speedytree] NJ trees functions uses
    // a N x N matrix. So to create this matrix, I generate all the combinations
    // using a cartesian product and then fill the matrix.
    let it = (0..n).cartesian_product(0..n);
    for (i, j) in it {
        match map.get(&(unique_names[i].clone(), unique_names[j].clone())) {
            Some(&dist) => {
                matrix[i][j] = dist;
            }
            None => {
                let dist = map
                    .get(&(unique_names[j].clone(), unique_names[i].clone()))
                    .unwrap();
                matrix[i][j] = *dist;
            }
        }
    }
    DistanceMatrix {
        matrix,
        names: unique_names,
    }
}

/// Write a PHYLIP file from a distance matrice
pub fn to_phylip(dist: DistanceMatrix, output: &str) -> anyhow::Result<()> {
    let mut path = PathBuf::from(output);
    path.push("distance.phylip");

    let mut file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;

    let mut mstr = String::with_capacity(2);
    mstr.push(char::from_digit(dist.names.len().try_into().unwrap(), 10).unwrap());
    mstr.push('\n');

    file.write_all(mstr.as_bytes())?;

    for i in 0..dist.names.len() {
        // Using format! here despite possible performance issue to limit decimal printed
        let val = dist.matrix[i].iter().format(" ").to_string();

        let mut mstr = String::with_capacity(dist.names[i].len() + val.len() + 2);
        mstr.push_str(&dist.names[i]);
        mstr.push(' ');
        mstr.push_str(&val);
        mstr.push('\n');

        file.write_all(mstr.as_bytes())?;
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
