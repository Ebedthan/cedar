// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::collections::HashMap;

use crate::dist::rescue::RESCUE_PVALUE_THRESHOLD;
use crate::mash::uncertainty::{mash_significance, DistanceEstimate, Reliability};
use finch::serialization::Sketch;

/// Build connectivity-analysis edges from a set of distance estimates.
/// Classifying a pairs as divergent if it fails EITHER the CI-based
/// reliability bar OR Mash's own significance test.
/// It is a "both checks" standard (attempt_kmer_rescue): a candidate k can't
/// pass connectivity purely by making noise look precise.
pub(crate) fn connectivity_edges(
    estimates: &[DistanceEstimate],
    sketches: &[Sketch],
    kmer_length: u8,
) -> Vec<crate::build::connectivity::Edge> {
    let sizes: HashMap<String, u64> = sketches
        .iter()
        .map(|s| (s.name.clone(), s.seq_length))
        .collect();

    estimates
        .iter()
        .map(|e| {
            let ci_divergent = matches!(
                e.reliability,
                Reliability::Unreliable | Reliability::NoSharedHashes
            );
            let query_size = *sizes.get(&e.query_path).unwrap_or(&0);
            let reference_size = *sizes.get(&e.reference_path).unwrap_or(&0);
            let pvalue = mash_significance(
                e.shared_hashes,
                e.total_hashes,
                query_size,
                reference_size,
                kmer_length,
            );
            let is_divergent = ci_divergent || pvalue > RESCUE_PVALUE_THRESHOLD;

            crate::build::connectivity::Edge {
                query: e.query.clone(),
                reference: e.reference.clone(),
                is_divergent,
                weight: e.relative_uncertainty.unwrap_or(0.0),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mash::distance::compute_distances;
    use crate::mash::uncertainty::mash_pvalue;
    use finch::serialization::SketchDistance;
    use std::io::Read;
    use tempfile;

    use itertools::Itertools;

    use std::fs;

    use crate::build::matrix::{distance_to_matrix, to_phylip};

    use speedytree::DistanceMatrix;

    use crate::dist::rescue::rescue_kmer_candidates;
    use crate::mash::uncertainty::{annotate_with_uncertainty, jaccard_ci_95};

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

    fn make_sketch_distance(
        query: &str,
        reference: &str,
        jaccard: f64,
        common_hashes: u64,
        total_hashes: u64,
    ) -> SketchDistance {
        let mash_distance = (-((2.0 * jaccard) / (1.0 + jaccard)).ln() / 21.0)
            .max(0.0)
            .min(1.0);
        SketchDistance {
            containment: jaccard, // not used by these tests
            jaccard,
            mash_distance,
            common_hashes,
            total_hashes,
            query: query.to_string(),
            reference: reference.to_string(),
        }
    }

    #[test]
    fn test_reliability_classification_matches_expected_regimes() {
        // Near-identical strains: tight Jaccard estimate -> reliable, even
        // though the resulting distance is tiny.
        let close = annotate_with_uncertainty(&make_sketch_distance("A", "B", 0.95, 950, 1000), 21);
        assert_eq!(close.reliability, Reliability::Reliable);
        assert!(close.relative_uncertainty.unwrap() < 0.10);

        // Same genus, different species: moderate relative uncertainty.
        let moderate =
            annotate_with_uncertainty(&make_sketch_distance("A", "C", 0.10, 100, 1000), 21);
        assert_eq!(moderate.reliability, Reliability::Borderline);

        // Cross-genus: very few shared hashes, estimate is mostly noise.
        let distant =
            annotate_with_uncertainty(&make_sketch_distance("A", "D", 0.005, 5, 1000), 21);
        assert_eq!(distant.reliability, Reliability::Unreliable);

        // No shared hashes at all: nothing to estimate uncertainty from.
        let none_shared =
            annotate_with_uncertainty(&make_sketch_distance("A", "E", 0.0, 0, 1000), 21);
        assert_eq!(none_shared.reliability, Reliability::NoSharedHashes);
        assert!(none_shared.relative_uncertainty.is_none());
    }

    #[test]
    fn test_distance_ci_is_within_bounds() {
        let estimate =
            annotate_with_uncertainty(&make_sketch_distance("A", "C", 0.10, 100, 1000), 21);
        assert!(estimate.mash_distance_ci_95_low <= estimate.mash_distance);
        assert!(estimate.mash_distance_ci_95_high >= estimate.mash_distance);
        assert!((0.0..=1.0).contains(&estimate.mash_distance_ci_95_low));
        assert!((0.0..=1.0).contains(&estimate.mash_distance_ci_95_high));

        // The exact CI on J should likewise bracket the point estimate.
        assert!(estimate.jaccard_ci_95_low <= estimate.jaccard);
        assert!(estimate.jaccard_ci_95_high >= estimate.jaccard);
    }

    #[test]
    fn test_clopper_pearson_boundary_cases() {
        // x=0: exact lower bound is 0, no Beta distribution needed for it.
        let (low, _high) = jaccard_ci_95(0, 1000).unwrap();
        assert_eq!(low, 0.0);

        // x=n: exact upper bound is 1.
        let (_low, high) = jaccard_ci_95(1000, 1000).unwrap();
        assert_eq!(high, 1.0);

        // total == 0: nothing to compute.
        assert!(jaccard_ci_95(0, 0).is_none());
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

    #[test]
    fn test_mash_significance_matches_paper_example() {
        // Ondov et al. state that s=400, true J=0.1 gives P(30<=x<=50) > 0.9.
        let p_ge_30 = 1.0 - mash_pvalue(30, 400, 0.1);
        let p_ge_51 = 1.0 - mash_pvalue(51, 400, 0.1);
        // p_ge_30 - p_ge_51 == P(30<=X<=50) is awkward via our p-value framing
        // (mash_pvalue(x,..) = P(X>=x)), so compute directly instead:
        let p_between = mash_pvalue(30, 400, 0.1) - mash_pvalue(51, 400, 0.1);
        assert!(p_between > 0.9);
    }

    #[test]
    fn test_rescue_kmer_candidates_always_includes_floor() {
        assert_eq!(rescue_kmer_candidates(21, 14), vec![19, 17, 15, 14]);
        assert_eq!(rescue_kmer_candidates(15, 14), vec![14]);
        assert_eq!(rescue_kmer_candidates(21, 21), Vec::<u8>::new());
        assert_eq!(rescue_kmer_candidates(14, 15), Vec::<u8>::new());
    }
}
