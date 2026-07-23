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

use statrs::distribution::{Beta, ContinuousCDF};

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

/// How much a pairwise distance estimate can be trusted, based on the
/// relative half-width of the 95% CI on the underlying Jaccard estimate
/// (not on the Mash distance itself, see the comment on
/// `compute_distances_with_uncertainty` for why).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reliability {
    /// The two sketches shared no hashes at all; there is no signal to put
    /// an uncertainty estimate on in the first place.
    NoSharedHashes,
    Reliable,
    Borderline,
    Unreliable,
}

impl Reliability {
    /// Relative CI half-width at/under which a pair is considered reliable.
    const RELIABLE_THRESHOLD: f64 = 0.10;

    /// Relative CI half-width at/under which a pair is considered borderline
    /// (above this, it's flagged unreliable).
    const BORDERLINE_THRESHOLD: f64 = 0.30;

    fn classify(relative_uncertainty: Option<f64>, common_hashes: u64) -> Self {
        if common_hashes == 0 {
            return Reliability::NoSharedHashes;
        }
        match relative_uncertainty {
            // No relative uncertainty could be computed (e.g. J == 0, or an
            // exact match with J == 1 where the CI collapses to a point),
            // nothing here to flag as unreliable.
            None => Reliability::Reliable,
            Some(r) if r <= Self::RELIABLE_THRESHOLD => Reliability::Reliable,
            Some(r) if r <= Self::BORDERLINE_THRESHOLD => Reliability::Borderline,
            _ => Reliability::Unreliable,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Reliability::NoSharedHashes => "no_shared_hashes",
            Reliability::Reliable => "reliable",
            Reliability::Borderline => "borderline",
            Reliability::Unreliable => "unreliable",
        }
    }
}

impl std::fmt::Display for Reliability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// A pairwise Mash distance annotated with an exact confidence interval,
/// as produced by `cedar dist` (no tree is built from these).
#[derive(Debug, Clone)]
pub struct DistanceEstimate {
    pub query: String,
    pub reference: String,
    pub jaccard: f64,
    pub jaccard_ci_95_low: f64,
    pub jaccard_ci_95_high: f64,
    pub mash_distance: f64,
    pub mash_distance_ci_95_low: f64,
    pub mash_distance_ci_95_high: f64,
    pub shared_hashes: u64,
    pub total_hashes: u64,

    /// Relative half-width of the 95% CI on the *Jaccard* estimate, the
    /// quantity `Reliability` is actually classified on. Reported here as
    /// a diagnostic, not as the primary "how far apart are these genomes"
    /// figure (that's `mash_distance` / `ci_95_*`).
    pub relative_uncertainty: Option<f64>,

    pub reliability: Reliability,
}

// Exact 95% confidence interval (Clopper-Pearson) for the Jaccard estimates (Ondov et al., 2016)
// Treating the shared-hash count as Binomial(s, J) following Ondov et al.,
// hypergeometric-to-binomial approximation Mash itself uses.
//
// We use the exact binomial CDF as Ondov et al., (Figure S1)
// `shared` is the observed shared-hash count (x), `total` is the effective comparison
// size (s, finch's `total_hashes`). Returns `None` if `total == 0`.
fn jaccard_ci_95(shared: u64, total: u64) -> Option<(f64, f64)> {
    if total == 0 {
        return None;
    }

    let x = shared as f64;
    let n = total as f64;
    let alpha = 0.05;

    // Clopper-Pearson: P(X >= x | n, p) = I_p(x, n-x+1), the regularized
    // incomplete beta function, so its bounds are inverse-beta quantiles.
    // Boundary cases (x = 0, x = n) are handled separately since the corresponding Beta
    // shape parameter would otherwise be 0 (invalid).
    let low = if shared == 0 {
        0.0
    } else {
        Beta::new(x, n - x + 1.0)
            .map(|b| b.inverse_cdf(alpha / 2.0))
            .unwrap_or(0.0)
    };

    let high = if shared == total {
        1.0
    } else {
        Beta::new(x + 1.0, n - x)
            .map(|b| b.inverse_cdf(1.0 - alpha / 2.0))
            .unwrap_or(1.0)
    };

    Some((low, high))
}

/// Mash distance as a function of the Jaccard index
/// Matching finch's `distance()` computation exactly
/// D = -(1 / k) * ln(2J / (1 + J)), clamped to [0, 1]
/// Used to transform the exact CI bounds on J directly
/// into bounds on D, rather than linearizing with the delta method
fn mash_distance_from_jaccard(jaccard: f64, kmer_length: u8) -> f64 {
    if jaccard <= 0.0 {
        return 1.0;
    }
    let k = kmer_length as f64;
    let d = -((2.0 * jaccard) / (1.0 + jaccard)).ln() / k;
    d.clamp(0.0, 1.0)
}

fn annotate_with_uncertainty(d: SketchDistance, kmer_length: u8) -> DistanceEstimate {
    let jaccard_ci = jaccard_ci_95(d.common_hashes, d.total_hashes);
    let (jaccard_ci_95_low, jaccard_ci_95_high) = jaccard_ci.unwrap_or((d.jaccard, d.jaccard));

    // D is monotonically *decreasing* in J, so the distance's upper bound
    // comes from J's lower bound, and vice versa.
    let (mash_distance_ci_95_low, mash_distance_ci_95_high) = match jaccard_ci {
        Some((j_low, j_high)) => (
            mash_distance_from_jaccard(j_high, kmer_length),
            mash_distance_from_jaccard(j_low, kmer_length),
        ),
        None => (d.mash_distance, d.mash_distance),
    };

    // Relative half-width of the exact CI on J
    let relative_uncertainty = jaccard_ci
        .filter(|_| d.jaccard > 0.0)
        .map(|(low, high)| (high - low) / (2.0 * d.jaccard));
    let reliability = Reliability::classify(relative_uncertainty, d.common_hashes);

    DistanceEstimate {
        query: extract_basename(&d.query),
        reference: extract_basename(&d.reference),
        jaccard: d.jaccard,
        jaccard_ci_95_low,
        jaccard_ci_95_high,
        mash_distance: d.mash_distance,
        mash_distance_ci_95_low,
        mash_distance_ci_95_high,
        shared_hashes: d.common_hashes,
        total_hashes: d.total_hashes,
        relative_uncertainty,
        reliability,
    }
}

/// Compute pairwise Mash distances annotated with an uncertainty estimate
/// on each one (the `cedar dist` command). Self-comparisons are dropped,
/// unlike `compute_distances`, which keeps them for `distance_to_matrix`'s
/// diagonal, a self-pair's distance is trivially 0 and not informative in
/// a distance-with-uncertainty table.
///
/// All sketches in a single cedar run share the same k-mer size (they're
/// generated together from the same CLI invocation), so k is read once
/// from the first sketch rather than tracked per pair.
pub fn compute_distances_with_uncertainty(sketches: Vec<Sketch>) -> Vec<DistanceEstimate> {
    let kmer_length = sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21);
    let distances = compute_distances(sketches);

    distances
        .into_iter()
        .filter(|d| extract_basename(&d.query) != extract_basename(&d.reference))
        .map(|d| annotate_with_uncertainty(d, kmer_length))
        .collect()
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

/// Eq. 1 (Ondov et al. 2016): probability that a specific k-mer appears
/// somewhere in a random genome of length `n`, alphabet size 4.
fn kmer_hit_probability(genome_size: u64, kmer_length: u8) -> f64 {
    1.0 - (1.0 - 4f64.powi(-(kmer_length as i32))).powf(genome_size as f64)
}

/// Eq. 5: expected Jaccard index between two *unrelated* random genomes
/// given each genome's own chance-hit probability from Eq. 1.
fn background_jaccard(p_x: f64, p_y: f64) -> f64 {
    let denom = p_x + p_y - p_x * p_y;
    if denom <= 0.0 {
        0.0
    } else {
        (p_x * p_y) / denom
    }
}

/// Eq. 8: P(shared hashes >= x | Binomial(s, r)) via the same
/// incomplete-beta identity used for the Clopper-Pearson CI.
fn mash_pvalue(x: u64, s: u64, r: f64) -> f64 {
    if x == 0 {
        return 1.0;
    }

    if r <= 0.0 {
        return 0.0;
    }

    if r >= 1.0 {
        return 1.0;
    }

    Beta::new(x as f64, (s - x + 1) as f64)
        .map(|b| b.cdf(r))
        .unwrap_or(1.0)
}

/// Mash's significance test:
/// Is the observed sharing between two genomes of sizes `size_x` / `size_y` distinguishable from chance?
pub fn mash_significance(
    shared: u64,
    total: u64,
    size_x: u64,
    size_y: u64,
    kmer_length: u8,
) -> f64 {
    let r = background_jaccard(
        kmer_hit_probability(size_x, kmer_length),
        kmer_hit_probability(size_y, kmer_length),
    );
    mash_pvalue(shared, total, r)
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
        let close = annotate_with_uncertainty(make_sketch_distance("A", "B", 0.95, 950, 1000), 21);
        assert_eq!(close.reliability, Reliability::Reliable);
        assert!(close.relative_uncertainty.unwrap() < 0.10);

        // Same genus, different species: moderate relative uncertainty.
        let moderate =
            annotate_with_uncertainty(make_sketch_distance("A", "C", 0.10, 100, 1000), 21);
        assert_eq!(moderate.reliability, Reliability::Borderline);

        // Cross-genus: very few shared hashes, estimate is mostly noise.
        let distant = annotate_with_uncertainty(make_sketch_distance("A", "D", 0.005, 5, 1000), 21);
        assert_eq!(distant.reliability, Reliability::Unreliable);

        // No shared hashes at all: nothing to estimate uncertainty from.
        let none_shared =
            annotate_with_uncertainty(make_sketch_distance("A", "E", 0.0, 0, 1000), 21);
        assert_eq!(none_shared.reliability, Reliability::NoSharedHashes);
        assert!(none_shared.relative_uncertainty.is_none());
    }

    #[test]
    fn test_distance_ci_is_within_bounds() {
        let estimate =
            annotate_with_uncertainty(make_sketch_distance("A", "C", 0.10, 100, 1000), 21);
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
}
