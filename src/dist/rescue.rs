use crate::mash::sketch::{create_sketches, k_computing};
use crate::mash::uncertainty::{
    annotate_with_uncertainty, compute_base_distances_with_uncertainty, jaccard_ci,
    mash_significance, DistanceEstimate, Reliability,
};
use finch::distance::distance;
use finch::serialization::Sketch;
use std::collections::HashMap;
use std::path::PathBuf;

/// Ceiling used only for the "would more sketching alone fix this" gate.
const RESCUE_SKETCH_SIZE_CEILING: u64 = 100_000;

/// k-mer size step for the rescue descent. Not yet exposed via CLI.
const RESCUE_KMER_STEP: u8 = 2;

/// Target probability for Fofanov's chance-collision formula,
/// used only to compute the rescue floor.
/// This is NOT the same thing as `--rescue-pvalue`: this one bounds
/// the chance of *any* random k-mer match given a genome size
/// (a k-selection question); `rescue_pvalue` bounds the chance of observing
/// *this many* shared hashes given a background rate (Mash's significance test).
/// They happen to share a conventional default of 0.01, but they are different
/// statistical objects and are not coupled to the same user-facing flag.
const FOFANOV_TARGET_PROBABILITY: f64 = 0.01;

/// Summary counts from a `compute_distances_with_uncertainty` run, used to
/// print the always-on stderr summary at the CLI boundary.
#[derive(Debug, Clone, Default)]
pub struct RescueSummary {
    pub needs_larger_sketch: usize,
    pub unreliable_rescued: usize,
    pub borderline_rescued: usize,
    pub still_unresolvable: usize,
}

/// Build the descending list of candidate k-mer sizes to try during a
/// rescue, from `current_kmer` down to and including `floor`. Always
/// includes the floor itself, even when the step size doesn't evenly
/// divide the gap, otherwise a candidate could be silently skipped.
pub(crate) fn rescue_kmer_candidates(current_kmer: u8, floor: u8) -> Vec<u8> {
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
/// below `target` (Borderline pairs get `target = RELIABLE_THRESHOLD`;
/// Unreliable/NoSharedHashes pairs get `target = BORDERLINE_THRESHOLD`),
/// by trying progressively smaller k-mer sizes.
/// Re-sketches only the two genomes in this pair (not the whole dataset)
/// at each candidate k, using the same fixed sketch size and seed as the
/// rest of the run. An existing sketch can't be "re-k-mer-ized" after the
/// fact: k is fixed at sketch-creation time.
///
/// A candidate k is accepted only if it passes BOTH:
///   - relative CI half-width on J at or under `target`, and
///   - Mash's own significance test at or under `rescue_pvalue`: the
///     observed sharing must be distinguishable from two random,
///     unrelated genomes at this k, not just numerically less noisy.
///
/// The floor is `k_computing` applied to the LARGER of the two genome
/// sizes: the larger genome is the more permissive one for chance k-mer
/// collisions, so it sets the k below which "success" would just be
/// k-mer-space exhaustion, not real homology.
#[allow(clippy::too_many_arguments)]
fn attempt_kmer_rescue(
    query_path: &str,
    reference_path: &str,
    query_size: u64,
    reference_size: u64,
    current_kmer: u8,
    sketch_size: usize,
    seed: u64,
    confidence: f64,
    rescue_pvalue: f64,
    target: f64,
) -> anyhow::Result<Option<DistanceEstimate>> {
    let floor = k_computing(
        query_size.max(reference_size) as u32,
        FOFANOV_TARGET_PROBABILITY,
    );
    let candidates = rescue_kmer_candidates(current_kmer, floor);

    let tmp = tempfile::tempdir()?;
    let tmp_path = tmp.path().to_string_lossy().to_string();
    let inputs = [PathBuf::from(query_path), PathBuf::from(reference_path)];

    for candidate_k in candidates {
        let sketch_paths = create_sketches(&inputs, candidate_k, sketch_size, seed, &tmp_path)?;

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

        let Some((low, high)) = jaccard_ci(d.common_hashes, d.total_hashes, confidence) else {
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

        if relative_uncertainty <= target && pvalue <= rescue_pvalue {
            let mut estimate = annotate_with_uncertainty(&d, candidate_k, confidence);
            estimate.rescued = true;
            return Ok(Some(estimate));
        }
    }

    Ok(None)
}

/// Search for the smallest sketch size in `[min_size, max_size]` whose
/// exact CI relative half-width, for a hypothetical pair with Jaccard
/// estimate `j`, would be at or under `target_precision`.
///
/// Uses the exact CI (`jaccard_ci`) at each candidate size rather than a
/// closed-form approximation. Relative uncertainty is monotonically
/// non-increasing in sketch size for fixed `j`, so binary search is valid.
///
/// Returns `None` if even `max_size` can't achieve the target — the gate
/// `compute_distances_with_uncertainty` uses to tell "needs more
/// sketching" apart from "genuinely needs a smaller k."
fn min_sketch_size_for_precision(
    j: f64,
    target_precision: f64,
    min_size: u64,
    max_size: u64,
    confidence: f64,
) -> Option<u64> {
    if j <= 0.0 {
        return None;
    }

    let meets_target = |s: u64| -> bool {
        let x = (j * s as f64).round() as u64;
        match jaccard_ci(x, s, confidence) {
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

/// Compute pairwise Mash distances annotated with an uncertainty estimate.
///
/// When `no_rescue` is false (the default), attempts to upgrade any pair
/// below its tier's target: Borderline pairs are pushed toward Reliable,
/// Unreliable/NoSharedHashes pairs are pushed toward Borderline-or-better.
/// Before spending a k-mer rescue attempt, checks whether a larger sketch
/// alone (up to a ceiling) would already clear the bar. Sketch size
/// governs sampling precision, not whether real signal exists, so growing
/// it is left to the user (`-s`), not done silently per pair.
///
/// When `no_rescue` is true, every pair is reported exactly as computed at
/// the run's single k-mer size, with no attempt to improve it.
pub fn compute_distances_with_uncertainty(
    sketches: Vec<Sketch>,
    sketch_size: usize,
    seed: u64,
    confidence: f64,
    rescue_pvalue: f64,
    no_rescue: bool,
) -> (Vec<DistanceEstimate>, RescueSummary) {
    // Captured before `sketches` is moved into compute_base_distances_with_uncertainty.
    let sizes: HashMap<String, u64> = sketches
        .iter()
        .map(|s| (s.name.clone(), s.seq_length))
        .collect();
    let kmer_length = sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21);

    let base_estimates = compute_base_distances_with_uncertainty(sketches, confidence);

    if no_rescue {
        return (base_estimates, RescueSummary::default());
    }

    let mut needs_larger_sketch = 0usize;
    let mut unreliable_rescued = 0usize;
    let mut borderline_rescued = 0usize;
    let mut still_unresolvable = 0usize;

    let estimates: Vec<DistanceEstimate> = base_estimates
        .into_iter()
        .map(|estimate| {
            let target = match estimate.reliability {
                Reliability::Reliable => return estimate,
                Reliability::Borderline => Reliability::RELIABLE_THRESHOLD,
                Reliability::Unreliable | Reliability::NoSharedHashes => {
                    Reliability::BORDERLINE_THRESHOLD
                }
            };
            let started_borderline = matches!(estimate.reliability, Reliability::Borderline);

            let fixable_by_sketch_size = min_sketch_size_for_precision(
                estimate.jaccard,
                target,
                sketch_size as u64,
                RESCUE_SKETCH_SIZE_CEILING,
                confidence,
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
                confidence,
                rescue_pvalue,
                target,
            ) {
                Ok(Some(rescued)) => {
                    if started_borderline {
                        borderline_rescued += 1;
                    } else {
                        unreliable_rescued += 1;
                    }
                    rescued
                }
                _ => {
                    still_unresolvable += 1;
                    estimate
                }
            }
        })
        .collect();

    (
        estimates,
        RescueSummary {
            needs_larger_sketch,
            unreliable_rescued,
            borderline_rescued,
            still_unresolvable,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use finch::filtering::FilterParams;
    use finch::sketch_schemes::{KmerCount, SketchParams};

    #[test]
    fn test_rescue_kmer_candidates_always_includes_floor() {
        assert_eq!(rescue_kmer_candidates(21, 14), vec![19, 17, 15, 14]);
        assert_eq!(rescue_kmer_candidates(15, 14), vec![14]);
        assert_eq!(rescue_kmer_candidates(21, 21), Vec::<u8>::new());
        assert_eq!(rescue_kmer_candidates(14, 15), Vec::<u8>::new());
    }

    #[test]
    fn test_min_sketch_size_for_precision() {
        let j = 0.10;
        let target = 0.10;
        let size = min_sketch_size_for_precision(j, target, 200, 100_000, 0.95).unwrap();
        let x = (j * size as f64).round() as u64;
        let (low, high) = jaccard_ci(x, size, 0.95).unwrap();
        let achieved = (high - low) / (2.0 * j);
        assert!(achieved <= target + 1e-6);

        let loose = min_sketch_size_for_precision(0.10, 0.30, 200, 100_000, 0.95).unwrap();
        let tight = min_sketch_size_for_precision(0.10, 0.05, 200, 100_000, 0.95).unwrap();
        assert!(tight >= loose);

        assert!(min_sketch_size_for_precision(0.0005, 0.10, 200, 100_000, 0.95).is_none());
    }

    fn make_sketch(
        name: &str,
        shared_range: std::ops::Range<u64>,
        unique_offset: u64,
        unique_count: u64,
    ) -> Sketch {
        let mut hashes: Vec<KmerCount> = shared_range
            .map(|h| KmerCount {
                hash: h,
                kmer: vec![],
                count: 1,
                extra_count: 0,
                label: None,
            })
            .collect();
        hashes.extend((0..unique_count).map(|i| KmerCount {
            hash: unique_offset + i,
            kmer: vec![],
            count: 1,
            extra_count: 0,
            label: None,
        }));
        hashes.sort_by_key(|k| k.hash);
        Sketch {
            name: name.to_string(),
            seq_length: 5_000_000,
            num_valid_kmers: hashes.len() as u64,
            comment: String::new(),
            hashes,
            filter_params: FilterParams::default(),
            sketch_params: SketchParams::Mash {
                kmers_to_sketch: 1000,
                final_size: 1000,
                no_strict: false,
                kmer_length: 21,
                hash_seed: 42,
            },
        }
    }

    #[test]
    fn test_pair_fixable_by_larger_sketch_never_attempts_kmer_rescue() {
        let p = make_sketch("P", 0..20, 1_000_000, 180);
        let q = make_sketch("Q", 0..20, 2_000_000, 180);

        let (estimates, summary) =
            compute_distances_with_uncertainty(vec![p, q], 200, 42, 0.95, 0.01, false);

        assert_eq!(estimates.len(), 1);
        assert_eq!(estimates[0].reliability, Reliability::Unreliable);
        assert!(!estimates[0].rescued);
        assert_eq!(summary.needs_larger_sketch, 1);
        assert_eq!(summary.unreliable_rescued, 0);
        assert_eq!(summary.borderline_rescued, 0);
    }

    #[test]
    fn test_genuinely_unresolvable_pair_fails_gracefully_without_panicking() {
        let r = make_sketch("R", 0..1, 10_000_000, 9999);
        let s = make_sketch("S", 0..1, 20_000_000, 9999);

        let (estimates, summary) =
            compute_distances_with_uncertainty(vec![r, s], 10000, 42, 0.95, 0.01, false);

        assert_eq!(estimates.len(), 1);
        assert_eq!(estimates[0].reliability, Reliability::Unreliable);
        assert!(!estimates[0].rescued);
        assert_eq!(summary.still_unresolvable, 1);
    }

    #[test]
    fn test_no_rescue_flag_skips_everything() {
        let r = make_sketch("R", 0..1, 10_000_000, 9999);
        let s = make_sketch("S", 0..1, 20_000_000, 9999);

        let (estimates, summary) =
            compute_distances_with_uncertainty(vec![r, s], 10000, 42, 0.95, 0.01, true);

        assert_eq!(estimates.len(), 1);
        assert_eq!(estimates[0].reliability, Reliability::Unreliable);
        // no_rescue must produce identical output to the base (unrescued)
        // computation -- summary is all zeros, no rescue attempted at all.
        assert_eq!(summary.needs_larger_sketch, 0);
        assert_eq!(summary.still_unresolvable, 0);
    }
}
