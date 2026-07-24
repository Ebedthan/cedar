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

use std::path::PathBuf;

/// k-mer size step for the rescue descent.
/// Not yet exposed via CLI.
const RESCUE_KMER_STEP: u8 = 2;

/// Mash's own significance threshold (Eq. 8)
/// Reused here as the acceptance bar for a rescue
const RESCUE_PVALUE_THRESHOLD: f64 = 0.01;

/// Ceiling used only for the "would more sketching alone fix this"
/// gate below
const RESCUE_SKETCH_SIZE_CEILING: u64 = 100_000;

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

    /// Raw scketch name preserved alongside the basename
    /// to ease rescue mechanism
    pub query_path: String,
    pub reference_path: String,

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

    /// The k-mer size these values were actually computed at.
    /// Equal to the run's k-mer size unless this pair was rescued
    pub kmer_size_used: u8,

    /// Whether this pair's values came from the rescue mechanism rather
    /// than the run's normal k-mer size
    pub rescued: bool,
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

/// Search for the smallest sketch size in `[min_size, max_size]` whose
/// exact Clopper-Pearson relative CI half-width, for a hypothetical pair
/// with Jaccard estimate `j`, would be at or under `target_precision`.
///
/// Uses the exact CI (`jaccard_ci_95`) at each candidate size rather than
/// a closed-form approximation. Relative uncertainty is monotonically
/// non-increasing in sketch size for fixed `j`, so binary search is valid.
///
/// Returns `None` if even `max_size` can't achieve the target. This is
/// the gate `compute_distances_with_uncertainty` uses to tell "needs more
/// sketching" apart from "genuinely needs a smaller k": sketch size
/// controls the *sampling precision* of the Jaccard estimate given some
/// true J, it cannot manufacture shared k-mers that aren't there.
fn min_sketch_size_for_precision(
    j: f64,
    target_precision: f64,
    min_size: u64,
    max_size: u64,
) -> Option<u64> {
    if j <= 0.0 {
        return None;
    }

    let meets_target = |s: u64| -> bool {
        let x = (j * s as f64).round() as u64;
        match jaccard_ci_95(x, s) {
            Some((low, high)) => (high - low) / (2.0 * j) <= target_precision,
            None => false,
        }
    };

    if !meets_target(max_size) {
        return None;
    }
    if meets_target(min_size) {
        return Some(min_size);
    }

    let (mut lo, mut hi) = (min_size, max_size);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if meets_target(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
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

fn annotate_with_uncertainty(d: &SketchDistance, kmer_length: u8) -> DistanceEstimate {
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
        query_path: d.query.clone(),
        reference_path: d.reference.clone(),
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
        kmer_size_used: kmer_length,
        rescued: false,
    }
}

/// Compute pairwise Mash distances annotated with an uncertainty
/// estimate, with NO rescue attempt. Used by `cedar build`'s
/// connectivity analysis, which needs every pair evaluated at a single,
/// uniform k-mer size, rescuing individual pairs to different k values
/// would break NJ's implicit assumption that every matrix entry is
/// comparable.
pub fn compute_base_distances_with_uncertainty(sketches: Vec<Sketch>) -> Vec<DistanceEstimate> {
    let kmer_length = sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21);
    compute_distances(sketches)
        .into_iter()
        .filter(|d| d.query != d.reference)
        .map(|d| annotate_with_uncertainty(&d, kmer_length))
        .collect()
}

/// Compute pairwise Mash distances annotated with an uncertainty estimate,
/// attempting a k-mer-size rescue (`attempt_kmer_rescue`) for any pair
/// that comes back Unreliable or NoSharedHashes at the run's k-mer size.
///
/// `sketch_size` and `seed` are needed here (not just `sketches`) because
/// rescuing a pair means re-sketching just its two genomes from their
/// original FASTA files at a smaller k.
pub fn compute_distances_with_uncertainty(
    sketches: Vec<Sketch>,
    sketch_size: usize,
    seed: u64,
) -> Vec<DistanceEstimate> {
    // SketchDistance doesn't carry genome sizes forward; keep them by
    // name for the rescue path's significance test and k-floor.
    let sizes: HashMap<String, u64> = sketches
        .iter()
        .map(|s| (s.name.clone(), s.seq_length))
        .collect();

    let kmer_length = sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21);
    let base_estimates = compute_base_distances_with_uncertainty(sketches);

    let mut needs_larger_sketch = 0usize;
    let mut rescued_count = 0usize;
    let mut still_unresolvable = 0usize;

    let estimates: Vec<DistanceEstimate> = base_estimates
        .into_iter()
        .map(|estimate| {
            if !matches!(
                estimate.reliability,
                Reliability::Unreliable | Reliability::NoSharedHashes
            ) {
                return estimate;
            }

            // Gate: would more sketching alone already fix this? If so,
            // don't spend a k-rescue attempt on it. Sketch size controls
            // sampling precision, not whether real signal exists, so this
            // is left to the user (`-s`) rather than auto-resketched.
            let fixable_by_sketch_size = min_sketch_size_for_precision(
                estimate.jaccard,
                Reliability::BORDERLINE_THRESHOLD,
                sketch_size as u64,
                RESCUE_SKETCH_SIZE_CEILING,
            )
            .is_some();

            if fixable_by_sketch_size {
                needs_larger_sketch += 1;
                return estimate;
            }

            let query_size = *sizes.get(&estimate.query_path).unwrap_or(&0);
            let reference_size = *sizes.get(&estimate.reference_path).unwrap_or(&0);

            match attempt_kmer_rescue(
                &estimate.query_path,
                &estimate.reference_path,
                query_size,
                reference_size,
                kmer_length,
                sketch_size,
                seed,
            ) {
                Ok(Some(rescued)) => {
                    rescued_count += 1;
                    rescued
                }
                _ => {
                    still_unresolvable += 1;
                    estimate
                }
            }
        })
        .collect();

    if needs_larger_sketch > 0 || rescued_count > 0 || still_unresolvable > 0 {
        println!(
            "Rescue summary: {} pair(s) resolvable with a larger sketch size (see relative_uncertainty), \
             {} pair(s) rescued with a smaller k-mer size (see kmer_size_used/rescued), \
             {} pair(s) remain unresolvable at any k down to the genome-size-appropriate floor.",
            needs_larger_sketch, rescued_count, still_unresolvable
        );
    }

    estimates
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

/// Build the descending list of candidate k-mer sizes to try during a
/// rescue, from `current_kmer` down to and including `floor`. Always
/// includes the floor itself, even when the step size doesn't evenly
/// divide the gap — otherwise a candidate could be silently skipped.
fn rescue_kmer_candidates(current_kmer: u8, floor: u8) -> Vec<u8> {
    let mut candidates = Vec::new();
    if current_kmer <= floor {
        return candidates;
    }
    let mut k = current_kmer.saturating_sub(RESCUE_KMER_STEP);
    while k > floor {
        candidates.push(k);
        k = k.saturating_sub(RESCUE_KMER_STEP);
    }
    candidates.push(floor);
    candidates
}

/// Attempt to rescue a pair whose reliability at the run's k-mer size is
/// Unreliable or NoSharedHashes, by trying progressively smaller k-mer
/// sizes. Re-sketches only the two genomes in this pair (not the whole
/// dataset) at each candidate k, using the same fixed sketch size and
/// seed as the rest of the run. An existing sketch can't be
/// "re-k-mer-ized" after the fact, k is fixed at sketch-creation time.
///
/// A candidate k is accepted only if it passes BOTH:
///   - relative CI half-width on J at or under the Borderline threshold
///     (the same acceptability bar used everywhere else), and
///   - Mash's own significance test: the observed sharing must be
///     distinguishable from two random, unrelated genomes at this k, not
///     just numerically less noisy.
///
/// The floor is `k_computing` applied to the LARGER of the two genome
/// sizes (at the existing 0.01 target probability): the larger genome is
/// the more permissive one for chance k-mer collisions, so it
/// sets the k below which "success" would just be k-mer-space
/// exhaustion, not real homology.
///
/// Returns `Ok(Some(estimate))` with `estimate.rescued = true` and
/// `estimate.kmer_size_used` set to the winning k, or `Ok(None)` if no
/// candidate down to the floor passes both checks.
fn attempt_kmer_rescue(
    query_path: &str,
    reference_path: &str,
    query_size: u64,
    reference_size: u64,
    current_kmer: u8,
    sketch_size: usize,
    seed: u64,
) -> anyhow::Result<Option<DistanceEstimate>> {
    let floor = super::sketch::k_computing(query_size.max(reference_size) as u32, 0.01);
    let candidates = rescue_kmer_candidates(current_kmer, floor);

    let tmp = tempfile::tempdir()?;
    let tmp_path = tmp.path().to_string_lossy().to_string();
    let inputs = [PathBuf::from(query_path), PathBuf::from(reference_path)];

    for candidate_k in candidates {
        let sketch_paths =
            super::sketch::create_sketches(&inputs, candidate_k, sketch_size, seed, &tmp_path)?;

        let mut candidate_sketches = Vec::new();
        for path in &sketch_paths {
            candidate_sketches.extend(finch::open_sketch_file(path)?);
        }
        if candidate_sketches.len() != 2 {
            continue;
        }

        let d = distance(&candidate_sketches[0], &candidate_sketches[1], false);
        let Ok(d) = d else {
            continue;
        };

        if d.jaccard <= 0.0 {
            continue;
        }

        let Some((low, high)) = jaccard_ci_95(d.common_hashes, d.total_hashes) else {
            continue;
        };
        let relative_uncertainty = (high - low) / (2.0 * d.jaccard);
        let pvalue = mash_significance(
            d.common_hashes,
            d.total_hashes,
            query_size,
            reference_size,
            candidate_k,
        );

        if relative_uncertainty <= Reliability::BORDERLINE_THRESHOLD
            && pvalue <= RESCUE_PVALUE_THRESHOLD
        {
            let mut estimate = annotate_with_uncertainty(&d, candidate_k);
            estimate.rescued = true;
            return Ok(Some(estimate));
        }
    }

    Ok(None)
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
