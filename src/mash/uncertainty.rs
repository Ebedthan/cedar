use crate::mash::distance::compute_distances;
use crate::mash::distance::extract_basename;
use finch::serialization::{Sketch, SketchDistance};
use statrs::distribution::{Beta, ContinuousCDF};

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
pub(crate) fn jaccard_ci_95(shared: u64, total: u64) -> Option<(f64, f64)> {
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

pub(crate) fn annotate_with_uncertainty(d: &SketchDistance, kmer_length: u8) -> DistanceEstimate {
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

/// Mash's significance test:
/// Is the observed sharing between two genomes of sizes `size_x` / `size_y` distinguishable from chance?
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
