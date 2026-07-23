pub mod build;
pub mod cli;
pub mod nwk;
pub mod utils;

use clap::Parser;

use anyhow::Result;

use crate::{build::build_tree_using_mash_distance, cli::Commands, utils::init_rayon_pool};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    init_rayon_pool(cli.threads);

    match cli.command {
        Commands::Build(args) => build_tree_using_mash_distance(&args, cli.threads),
    }
}
