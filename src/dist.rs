use std::collections::HashMap;
use std::path::Path;

use finch::{
    distance::distance,
    serialization::{Sketch, SketchDistance},
};
use itertools::Itertools;
use speedytree::DistanceMatrix;

// Function to compute distances
// Adapated from finch
// Permalink: https://github.com/onecodex/finch-rs/blob/47850b0cf6c506ef7b3f30966504f8732a1b888f/cli/src/main.rs#L315
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

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct Sequence {
    query: String,
    reference: String,
}

/// Transform a Vec of [finch]'s SketchDistance into a [speedytree] DistanceMatrix
pub fn sketches_distance_to_matrix(distances: Vec<SketchDistance>) -> DistanceMatrix {
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

    // Hold on! Here something is happening first the SketchDistance struct
    // DOES NOT CONTAINS all pairwise distances but only non repeating
    // distances. However, the [speedytree] NJ trees functions uses
    // a n x n matrix. So to create this matrix, I generate all the combinations
    // using a cartesian product and then filling the matrix.
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
