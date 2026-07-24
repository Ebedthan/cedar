use finch::distance::distance;
use finch::serialization::Sketch;
use finch::serialization::SketchDistance;
use itertools::Itertools;
use std::path::Path;

/// Cached basename extraction
pub fn extract_basename(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_string()
}

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
