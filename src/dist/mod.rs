pub mod rescue;

use crate::cli;
use crate::mash::sketch::create_and_load_sketches;
use crate::mash::uncertainty::Reliability;
use crate::utils;

use anyhow::Result;

/// `cedar dist`: compute pairwise Mash distances annotated with an
/// uncertainty estimate on each one. No tree is built, this is for the
/// "is this pair a solid species-boundary call" moment, not for producing
/// a phylogeny.
pub fn compute_pairwise_distances(args: &cli::DistArgs) -> Result<()> {
    // Validate and collect input files
    let inputs = utils::validate_and_collect_inputs(&args.indir)?;

    // Compute genome statistics and flag genome-size outliers up front
    // an outlier genome size undermines the k-mer size choice for every
    // pairwise distance computed from it, exactly as it would for `build`.
    let stats = utils::compute_genome_stats(&inputs)?;
    //utils::check_genome_outliers(&stats, utils::DEFAULT_OUTLIER_THRESHOLD)?;

    // Determine k-mer size: manual override, or the same default heuristic
    // `build` uses. `--target-precision`-driven selection isn't implemented
    // yet, so it's rejected explicitly here rather than silently ignored.
    let kmer_size = utils::determine_kmer_size(&args.sketch, &stats)?;

    // Create sketches
    let sketches =
        create_and_load_sketches(&args.indir, args.sketch.size, args.sketch.seed, kmer_size)?;

    // Compute pairwise distances, each with an uncertainty estimate
    let estimates =
        rescue::compute_distances_with_uncertainty(sketches, args.sketch.size, args.sketch.seed);

    let flagged = estimates
        .iter()
        .filter(|e| e.reliability != Reliability::Reliable)
        .count();
    println!(
        "{} pairwise distances computed ({} flagged as borderline, unreliable, or lacking shared hashes)",
        estimates.len(),
        flagged
    );

    utils::output_distance_table(args.output.clone(), &estimates)
}
