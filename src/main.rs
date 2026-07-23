// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

pub mod build;
pub mod cli;
pub mod nwk;
pub mod utils;

use clap::Parser;

use anyhow::Result;

use crate::{
    build::{build_tree_using_mash_distance, compute_pairwise_distances},
    cli::Commands,
    utils::init_rayon_pool,
};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    init_rayon_pool(cli.threads);

    match cli.command {
        Commands::Build(args) => build_tree_using_mash_distance(&args, cli.threads),
        Commands::Dist(args) => compute_pairwise_distances(&args),
    }
}
