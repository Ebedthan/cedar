use crate::mash::sketch::{create_sketches, k_computing};
use crate::mash::uncertainty::{
    annotate_with_uncertainty, compute_base_distances_with_uncertainty, jaccard_ci_95,
    mash_significance, DistanceEstimate, Reliability,
};
use finch::distance::distance;
use finch::serialization::Sketch;
use std::collections::HashMap;
use std::path::PathBuf;

/// Ceiling used only for the "would more sketching alone fix this"
/// gate below
const RESCUE_SKETCH_SIZE_CEILING: u64 = 100_000;

/// k-mer size step for the rescue descent.
/// Not yet exposed via CLI.
const RESCUE_KMER_STEP: u8 = 2;

/// Mash's own significance threshold (Eq. 8)
/// Reused here as the acceptance bar for a rescue
pub(crate) const RESCUE_PVALUE_THRESHOLD: f64 = 0.01;

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
    let floor = k_computing(query_size.max(reference_size) as u32, 0.01);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rescue_kmer_candidates_always_includes_floor() {
        assert_eq!(rescue_kmer_candidates(21, 14), vec![19, 17, 15, 14]);
        assert_eq!(rescue_kmer_candidates(15, 14), vec![14]);
        assert_eq!(rescue_kmer_candidates(21, 21), Vec::<u8>::new());
        assert_eq!(rescue_kmer_candidates(14, 15), Vec::<u8>::new());
    }
}
