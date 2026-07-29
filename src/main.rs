pub mod cli;
pub mod nwk;
pub mod utils;

pub mod build;
pub mod dist;
pub mod mash;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    utils::init_rayon_pool(cli.threads);

    match cli.command {
        cli::Commands::Build(args) => build::build_tree_using_mash_distance(&args, cli.threads),
        cli::Commands::Dist(args) => dist::compute_pairwise_distances(&args, cli.verbose),
        cli::Commands::Sketch(args) => {
            let inputs = utils::validate_and_collect_inputs(&args.inputs)?;
            let stats = utils::compute_genome_stats(&inputs)?;
            let kmer = utils::determine_kmer_size(&args.sketch, &stats)?;
            mash::sketch::sketch_only(
                &inputs,
                args.sketch.size,
                args.sketch.seed,
                kmer.kmer_size,
                &args.output_dir,
            )?;
            Ok(())
        }
    }
}
