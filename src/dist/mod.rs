use anyhow::Result;

use crate::cli;
use crate::mash::sketch::create_and_load_sketches;
use crate::mash::uncertainty::Reliability;
use crate::utils;

pub mod rescue;

/// `cedar dist`: compute pairwise Mash distances annotated with an
/// uncertainty estimate on each one.
/// No tree is built. This is for the  "is this pair a solid species-boundary call"
/// moment, not for producing a phylogeny. A summary is always printed to stderr.
/// The full TSV is only written if `-o/--output` is given.
pub fn compute_pairwise_distances(args: &cli::DistArgs) -> Result<()> {
    if !(args.uncertainty.confidence > 0.5 && args.uncertainty.confidence < 1.0) {
        anyhow::bail!(
            "--confidence must be strictly between 0.5 and 1.0, got {}",
            args.uncertainty.confidence
        );
    }
    if !(args.uncertainty.rescue_pvalue > 0.0 && args.uncertainty.rescue_pvalue < 1.0) {
        anyhow::bail!(
            "--rescue-pvalue must be strictly between 0.0 and 1.0, got {}",
            args.uncertainty.rescue_pvalue
        );
    }

    let inputs = utils::validate_and_collect_inputs(&args.inputs)?;

    let stats = utils::compute_genome_stats(&inputs)?;
    utils::check_genome_outliers(&stats, utils::DEFAULT_OUTLIER_THRESHOLD)?;

    let kmer_size = utils::determine_kmer_size(&args.sketch, &stats)?;

    let sketches =
        create_and_load_sketches(&inputs, args.sketch.size, args.sketch.seed, kmer_size)?;

    let (estimates, summary) = rescue::compute_distances_with_uncertainty(
        sketches,
        args.sketch.size,
        args.sketch.seed,
        args.uncertainty.confidence,
        args.uncertainty.rescue_pvalue,
        args.no_rescue,
    );

    let flagged = estimates
        .iter()
        .filter(|e| e.reliability != Reliability::Reliable)
        .count();

    eprintln!(
        "{} pairwise distances computed ({} flagged as borderline, unreliable, or lacking shared hashes)",
        estimates.len(),
        flagged
    );
    if !args.no_rescue {
        eprintln!(
            "Rescue: {} unreliable/no-shared-hashes => borderline-or-better, \
             {} borderline => reliable, {} pair(s) resolvable with a larger sketch size, \
             {} pair(s) still unresolvable at any k down to the genome-size-appropriate floor.",
            summary.unreliable_rescued,
            summary.borderline_rescued,
            summary.needs_larger_sketch,
            summary.still_unresolvable
        );
    }
    eprintln!(
        "Confidence level: {:.0}%, rescue significance threshold: {}",
        args.uncertainty.confidence * 100.0,
        args.uncertainty.rescue_pvalue
    );

    if let Some(output_path) = &args.output {
        utils::output_distance_table(Some(output_path.clone()), &estimates)?;
    }

    Ok(())
}
