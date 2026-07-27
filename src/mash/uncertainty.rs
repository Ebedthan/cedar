use crate::mash::distance::compute_distances;
use crate::mash::distance::extract_basename;

use finch::serialization::{Sketch, SketchDistance};
use statrs::distribution::{Beta, ContinuousCDF};

use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// How much a pairwise distance estimate can be trusted ?
/// Based on the relative half-width of the confidence interval on
/// the underlying Jaccard estimate (not on the Mash distance itself,
/// see the comment on `mash_distance_from_jaccard` for why). These
/// thresholds are fixed regardless of the confidence level requested.
/// They bucket the *relative* uncertainty, not the raw interval width,
/// so tightening or loosening `--confidence` doesn't change what
/// "reliable" means, only how hard it is to reach.
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
    pub(crate) const RELIABLE_THRESHOLD: f64 = 0.10;

    /// Relative CI half-width at/under which a pair is considered borderline
    /// (above this, it's flagged unreliable).
    pub const BORDERLINE_THRESHOLD: f64 = 0.30;

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

    /// Raw sketch name preserved alongside the basename, to ease the
    /// rescue mechanism's lookup of the original FASTA file.
    pub query_path: String,
    pub reference_path: String,

    pub jaccard: f64,
    pub jaccard_ci_low: f64,
    pub jaccard_ci_high: f64,
    pub mash_distance: f64,
    pub mash_distance_ci_low: f64,
    pub mash_distance_ci_high: f64,
    pub shared_hashes: u64,
    pub total_hashes: u64,

    /// Relative half-width of the CI on the *Jaccard* estimate, the
    /// quantity `Reliability` is actually classified on. Reported here as
    /// a diagnostic, not as the primary "how far apart are these genomes"
    /// figure (that's `mash_distance` / `mash_distance_ci_*`).
    pub relative_uncertainty: Option<f64>,

    pub reliability: Reliability,

    /// The k-mer size these values were actually computed at. Equal to
    /// the run's k-mer size unless this pair was rescued.
    pub kmer_size_used: u8,

    /// Whether this pair's values came from the rescue mechanism rather
    /// than the run's normal k-mer size.
    pub rescued: bool,
}

/// Exact confidence interval (Clopper-Pearson) for the Jaccard estimate
/// (Ondov et al., 2016), treating the shared-hash count as Binomial(s, J)
/// following the same hypergeometric-to-binomial approximation Mash itself
/// uses.
///
/// `shared` is the observed shared-hash count (x), `total` is the effective
/// comparison size (s, finch's `total_hashes`), and `confidence` is the
/// requested confidence level (e.g. 0.95 for a 95% interval, 0.99 for 99%).
/// Tightening `confidence` widens the interval for the same data.
/// Returns `None` if `total == 0`.
pub(crate) fn jaccard_ci(shared: u64, total: u64, confidence: f64) -> Option<(f64, f64)> {
    if total == 0 {
        return None;
    }

    let x = shared as f64;
    let n = total as f64;
    let alpha = 1.0 - confidence;

    // Clopper-Pearson: P(X >= x | n, p) = I_p(x, n - x + 1), the regularized
    // incomplete beta function, so its bounds are inverse-beta quantiles.
    // Boundary cases (x = 0, x = n) are handled separately since the
    // corresponding Beta shape parameter would otherwise be 0 (invalid).
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

pub(crate) fn annotate_with_uncertainty(
    d: &SketchDistance,
    kmer_length: u8,
    confidence: f64,
) -> DistanceEstimate {
    let jaccard_ci = jaccard_ci(d.common_hashes, d.total_hashes, confidence);
    let (jaccard_ci_low, jaccard_ci_high) = jaccard_ci.unwrap_or((d.jaccard, d.jaccard));

    // D is monotonically *decreasing* in J, so the distance's upper bound
    // comes from J's lower bound, and vice versa.
    let (mash_distance_ci_low, mash_distance_ci_high) = match jaccard_ci {
        Some((j_low, j_high)) => (
            mash_distance_from_jaccard(j_high, kmer_length),
            mash_distance_from_jaccard(j_low, kmer_length),
        ),
        None => (d.mash_distance, d.mash_distance),
    };

    // Relative half-width of the CI on J.
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
        jaccard_ci_low,
        jaccard_ci_high,
        mash_distance: d.mash_distance,
        mash_distance_ci_low,
        mash_distance_ci_high,
        shared_hashes: d.common_hashes,
        total_hashes: d.total_hashes,
        relative_uncertainty,
        reliability,
        kmer_size_used: kmer_length,
        rescued: false,
    }
}

/// Compute pairwise Mash distances annotated with an uncertainty
/// estimate, with NO rescue attempt. Used by `cedar build`'s connectivity
/// analysis, which needs every pair evaluated at a single, uniform k-mer
/// size, rescuing individual pairs to different k values would break
/// NJ's implicit assumption that every matrix entry is comparable.
pub fn compute_base_distances_with_uncertainty(
    sketches: Vec<Sketch>,
    confidence: f64,
) -> Vec<DistanceEstimate> {
    let kmer_length = sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21);
    compute_distances(sketches)
        .into_par_iter()
        .filter(|d| d.query != d.reference)
        .map(|d| annotate_with_uncertainty(&d, kmer_length, confidence))
        .collect()
}

/// Probability that a specific k-mer appears somewhere in a random
/// genome of length `n`, alphabet size 4.
/// Eq. 1 (Ondov et al. 2016)
fn kmer_hit_probability(genome_size: u64, kmer_length: u8) -> f64 {
    1.0 - (1.0 - 4f64.powi(-(kmer_length as i32))).powf(genome_size as f64)
}

/// Expected Jaccard index between two *unrelated* random genomes
/// given each genome's own chance-hit probability from `kmer_hit_probability`.
/// Eq. 5 (Ondov et al. 2016)
fn background_jaccard(p_x: f64, p_y: f64) -> f64 {
    let denom = p_x + p_y - p_x * p_y;
    if denom <= 0.0 {
        0.0
    } else {
        (p_x * p_y) / denom
    }
}

/// P(shared hashes >= x | Binomial(s, r)) via the same
/// incomplete-beta identity used for the Clopper-Pearson CI.
/// Eq. 8 (Ondov et al. 2016)
pub(crate) fn mash_pvalue(x: u64, s: u64, r: f64) -> f64 {
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

/// Mash's significance test: is the observed sharing between two genomes
/// of sizes `size_x` / `size_y` distinguishable from chance?
/// Returns a raw P-value; callers compare it against whatever threshold
/// they've been given (e.g. `--rescue-pvalue`), not a fixed constant here.
pub(crate) fn mash_significance(
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

/// Mash distance as a function of the Jaccard index, matching finch's
/// `distance()` computation exactly: D = -(1/k) * ln(2J / (1+J)), clamped
/// to [0, 1]. Used to transform the exact CI bounds on J directly into
/// bounds on D, rather than linearizing with the delta method.
fn mash_distance_from_jaccard(jaccard: f64, kmer_length: u8) -> f64 {
    if jaccard <= 0.0 {
        return 1.0;
    }
    let k = kmer_length as f64;
    let d = -((2.0 * jaccard) / (1.0 + jaccard)).ln() / k;
    d.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

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
            containment: jaccard,
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
        let close =
            annotate_with_uncertainty(&make_sketch_distance("A", "B", 0.95, 950, 1000), 21, 0.95);
        assert_eq!(close.reliability, Reliability::Reliable);
        assert!(close.relative_uncertainty.unwrap() < 0.10);

        let moderate =
            annotate_with_uncertainty(&make_sketch_distance("A", "C", 0.10, 100, 1000), 21, 0.95);
        assert_eq!(moderate.reliability, Reliability::Borderline);

        let distant =
            annotate_with_uncertainty(&make_sketch_distance("A", "D", 0.005, 5, 1000), 21, 0.95);
        assert_eq!(distant.reliability, Reliability::Unreliable);

        let none_shared =
            annotate_with_uncertainty(&make_sketch_distance("A", "E", 0.0, 0, 1000), 21, 0.95);
        assert_eq!(none_shared.reliability, Reliability::NoSharedHashes);
        assert!(none_shared.relative_uncertainty.is_none());
    }

    #[test]
    fn test_distance_ci_is_within_bounds() {
        let estimate =
            annotate_with_uncertainty(&make_sketch_distance("A", "C", 0.10, 100, 1000), 21, 0.95);
        assert!(estimate.mash_distance_ci_low <= estimate.mash_distance);
        assert!(estimate.mash_distance_ci_high >= estimate.mash_distance);
        assert!((0.0..=1.0).contains(&estimate.mash_distance_ci_low));
        assert!((0.0..=1.0).contains(&estimate.mash_distance_ci_high));

        assert!(estimate.jaccard_ci_low <= estimate.jaccard);
        assert!(estimate.jaccard_ci_high >= estimate.jaccard);
    }

    #[test]
    fn test_clopper_pearson_boundary_cases() {
        let (low, _high) = jaccard_ci(0, 1000, 0.95).unwrap();
        assert_eq!(low, 0.0);

        let (_low, high) = jaccard_ci(1000, 1000, 0.95).unwrap();
        assert_eq!(high, 1.0);

        assert!(jaccard_ci(0, 0, 0.95).is_none());
    }

    #[test]
    fn test_tighter_confidence_widens_the_interval() {
        // The whole point of --confidence: asking for 99% instead of 95%
        // must never produce a narrower interval for the same data.
        let (low95, high95) = jaccard_ci(100, 1000, 0.95).unwrap();
        let (low99, high99) = jaccard_ci(100, 1000, 0.99).unwrap();
        assert!(low99 <= low95);
        assert!(high99 >= high95);
    }

    #[test]
    fn test_mash_significance_matches_paper_example() {
        // Ondov et al. state that s=400, true J=0.1 gives P(30<=x<=50) > 0.9.
        let p_between = mash_pvalue(30, 400, 0.1) - mash_pvalue(51, 400, 0.1);
        assert!(p_between > 0.9);
    }

    #[test]
    fn test_kmer_hit_probability_and_background_jaccard_realistic_genomes() {
        // 5 Mbp genome at k=21: values already verified against a real
        // computation several turns back.
        let p = kmer_hit_probability(5_000_000, 21);
        assert!((p - 1.137e-6).abs() < 1e-9, "unexpected P(K): {}", p);

        // Two such genomes: background Jaccard should be roughly half of p
        // (since denom ~= 2p for tiny p).
        let r = background_jaccard(p, p);
        assert!(
            (r - 5.684e-7).abs() < 1e-10,
            "unexpected background r: {}",
            r
        );
    }

    #[test]
    fn test_kmer_hit_probability_increases_with_genome_size_and_decreases_with_k() {
        // Eq. 1's whole point: larger genomes and smaller k both make a
        // chance k-mer match more likely.
        assert!(kmer_hit_probability(10_000_000, 21) > kmer_hit_probability(5_000_000, 21));
        assert!(kmer_hit_probability(5_000_000, 15) > kmer_hit_probability(5_000_000, 21));
    }

    #[test]
    fn test_background_jaccard_symmetric_and_bounded() {
        let r1 = background_jaccard(0.001, 0.002);
        let r2 = background_jaccard(0.002, 0.001);
        assert_eq!(
            r1, r2,
            "background Jaccard must be symmetric in its two arguments"
        );
        assert!((0.0..=1.0).contains(&r1));

        // Degenerate case: both probabilities zero -> denom is zero -> defined as 0, not NaN.
        assert_eq!(background_jaccard(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_mash_distance_from_jaccard_matches_finch_formula() {
        // These exact values were confirmed against a real compiled
        // finch::distance() call earlier this conversation (A-vs-B and
        // A-vs-E pairs at k=21).
        assert!((mash_distance_from_jaccard(0.95, 21) - 0.0012).abs() < 1e-4);
        assert!((mash_distance_from_jaccard(0.005, 21) - 0.2195).abs() < 1e-4);

        // Identical sketches: distance is exactly 0.
        assert_eq!(mash_distance_from_jaccard(1.0, 21), 0.0);

        // No shared signal at all: distance is exactly 1 (fully diverged, by
        // convention), not an unclamped +infinity from ln(0).
        assert_eq!(mash_distance_from_jaccard(0.0, 21), 1.0);
        assert_eq!(mash_distance_from_jaccard(-0.5, 21), 1.0); // defensive: shouldn't occur, but must not panic
    }

    #[test]
    fn test_reliability_classify_exact_threshold_boundaries() {
        // Exactly at the boundary is inclusive on the tighter side (<=), not
        // exclusive -- worth locking in explicitly since the existing
        // regime tests (J=0.95/0.10/0.005) never land precisely on a cutoff.
        assert_eq!(
            Reliability::classify(Some(0.10), 100),
            Reliability::Reliable
        );
        assert_eq!(
            Reliability::classify(Some(0.100001), 100),
            Reliability::Borderline
        );
        assert_eq!(
            Reliability::classify(Some(0.30), 100),
            Reliability::Borderline
        );
        assert_eq!(
            Reliability::classify(Some(0.300001), 100),
            Reliability::Unreliable
        );

        // common_hashes == 0 always wins, regardless of relative_uncertainty.
        assert_eq!(
            Reliability::classify(Some(0.01), 0),
            Reliability::NoSharedHashes
        );
        assert_eq!(Reliability::classify(None, 0), Reliability::NoSharedHashes);
    }

    #[test]
    fn test_reliability_display_and_as_str() {
        assert_eq!(Reliability::Reliable.as_str(), "reliable");
        assert_eq!(format!("{}", Reliability::Borderline), "borderline");
        assert_eq!(format!("{}", Reliability::Unreliable), "unreliable");
        assert_eq!(
            format!("{}", Reliability::NoSharedHashes),
            "no_shared_hashes"
        );
    }

    #[test]
    fn test_mash_pvalue_edge_cases() {
        // No shared hashes at all: P-value is trivially 1 (can't reject the
        // null with zero observed signal).
        assert_eq!(mash_pvalue(0, 1000, 0.1), 1.0);

        // Background rate of exactly 0: any shared hash is maximally
        // significant.
        assert_eq!(mash_pvalue(5, 1000, 0.0), 0.0);

        // Background rate of 1 (degenerate): nothing is distinguishable from
        // chance.
        assert_eq!(mash_pvalue(5, 1000, 1.0), 1.0);
    }
}
