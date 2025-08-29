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
    build::{build_tree_from_genomes, build_tree_from_orthologous_groups},
    cli::Commands,
    utils::init_rayon_pool,
};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    init_rayon_pool(cli.threads);

    match cli.command {
        Commands::Build(args) => {
            if let Some(_) = args.genomes {
                build_tree_from_genomes(&args, cli.threads)
            } else {
                build_tree_from_orthologous_groups(&args, cli.threads)
            }
        }
        Commands::Compare(args) => compare_command(args),
    }
}

fn compare_command(args: cli::CompareArgs) -> anyhow::Result<()> {
    // TODO: implement functionality
    println!("Compare command not yet implemented: {:?}", args);
    Ok(())
}
