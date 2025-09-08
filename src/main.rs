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
    build::{build_tree_from_orthologous_groups, build_tree_using_mash_distance},
    cli::Commands,
    utils::init_rayon_pool,
};

fn main() -> Result<()> {
    let cli = cli::Cli::parse();
    init_rayon_pool(cli.threads);

    match cli.command {
        Commands::Build(args) => {
            if args.mash {
                build_tree_using_mash_distance(&args, cli.threads)?
            } else if args.ortholog {
                utils::check_required_tools()?;
                build_tree_from_orthologous_groups(&args, cli.threads)?
            }
            Ok(())
        }
        Commands::Compare(args) => compare_command(args),
    }
}

fn compare_command(args: cli::CompareArgs) -> anyhow::Result<()> {
    // TODO: implement functionality
    println!("Compare command not yet implemented: {:?}", args);
    Ok(())
}
