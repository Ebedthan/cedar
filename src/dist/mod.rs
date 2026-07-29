use anyhow::Result;

use crate::cli;
use crate::mash::sketch::{create_and_load_sketches, load_sketches};
use crate::utils;
use crate::utils::KmerSelection;

pub mod rescue;

/// `cedar dist`: compute pairwise Mash distances annotated with an
/// uncertainty estimate on each one.
/// No tree is built. This is for the  "is this pair a solid species-boundary call"
/// moment, not for producing a phylogeny. A summary is always printed to stderr.
/// The full TSV is only written if `-o/--output` is given.
pub fn compute_pairwise_distances(args: &cli::DistArgs, cli_verbose: bool) -> Result<()> {
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

    let inputs = if args.from_sketches.is_some() {
        Vec::new()
    } else {
        utils::validate_and_collect_inputs(&args.inputs)?
    };

    let (sketches, stats, kmer, outliers) = match &args.from_sketches {
        Some(sketch_dir) => {
            let sketches = load_sketches(sketch_dir)?;
            let stats: Vec<(String, usize)> = sketches
                .iter()
                .map(|s| (s.name.clone(), s.seq_length as usize))
                .collect();

            let kmer = KmerSelection {
                kmer_size: sketches.first().map(|s| s.sketch_params.k()).unwrap_or(21),
                user_specified: true, // true because loaded k as 'user-determined'
                mean_genome_size: None,
            };
            let outliers = utils::check_genome_outliers(&stats, utils::DEFAULT_OUTLIER_THRESHOLD);
            (sketches, stats, kmer, outliers)
        }
        None => {
            let stats = utils::compute_genome_stats(&inputs)?;

            let kmer = utils::determine_kmer_size(&args.sketch, &stats)?;
            let sketches = create_and_load_sketches(
                &inputs,
                args.sketch.size,
                args.sketch.seed,
                kmer.kmer_size,
            )?;
            let outliers = utils::check_genome_outliers(&stats, utils::DEFAULT_OUTLIER_THRESHOLD);
            (sketches, stats, kmer, outliers)
        }
    };

    let (estimates, summary) = rescue::compute_distances_with_uncertainty(
        sketches,
        args.sketch.size,
        args.sketch.seed,
        args.uncertainty.confidence,
        args.uncertainty.rescue_pvalue,
        args.no_rescue,
    );

    if let Some(output_path) = &args.output {
        utils::output_distance_table(
            Some(output_path.clone()),
            &estimates,
            args.uncertainty.confidence,
            args.uncertainty.rescue_pvalue,
        )?;
    }

    utils::print_run_summary(
        &stats,
        &outliers,
        &kmer,
        &args.sketch,
        &args.uncertainty,
        &estimates,
        if args.no_rescue { None } else { Some(&summary) },
        &args.output,
        cli_verbose,
    );

    Ok(())
}
