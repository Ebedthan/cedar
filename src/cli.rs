// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
    name = "cedar",
    about = "a phylogenomic toolkit",
    author,
    version,
    arg_required_else_help = true
)]
pub struct Cli {
    /// Number of threads to use
    #[arg(short, default_value_t = 1, value_name = "INT")]
    pub threads: usize,

    /// Add verbosity to program
    #[arg(short = 'v', long, action = clap::ArgAction::SetTrue)]
    pub verbose: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Phylogenetic tree building
    Build(BuildArgs),

    /// Compare trees
    Compare(CompareArgs),
}
#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Input directory containing FASTA files
    #[arg(value_name = "DIR")]
    pub indir: String,

    /// Build tree using mash distance
    #[arg(long, conflicts_with = "ortholog", help_heading = "Tree building mode")]
    pub mash: bool,

    /// Treat input files as orthologous groups
    #[arg(long, help_heading = "Tree building mode")]
    pub ortholog: bool,

    /// Output tree (Newick format) to FILE
    #[arg(short, value_name = "FILE")]
    pub output: Option<String>,

    /// Intermediate files directory.
    #[arg(long, value_name = "DIR")]
    pub tempdir: Option<String>,

    /// Bootstrap replicates
    #[arg(
        short = 'B',
        long = "boot",
        value_name = "INT",
        help_heading = "Tree options"
    )]
    pub bootstrap: Option<usize>,

    /// Jackknife replicates
    #[arg(
        short = 'J',
        long = "jack",
        value_name = "INT",
        help_heading = "Tree options"
    )]
    pub jacknife: Option<usize>,

    /// Subsampling proportion for jackknife
    #[arg(
        long,
        value_name = "FLOAT",
        default_value_t = 0.5,
        help_heading = "Tree options"
    )]
    pub jacknife_prop: f64,

    /// Compute canonical NJ tree
    #[arg(long, help_heading = "Tree options")]
    pub canonical: bool,

    /// Sketch size
    #[arg(
        short = 's',
        long,
        default_value_t = 1000,
        value_name = "INT",
        help_heading = "Sketching options"
    )]
    pub size: usize,

    /// Seed for hash function
    #[arg(
        short = 'S',
        long,
        default_value_t = 42,
        value_name = "INT",
        help_heading = "Sketching options"
    )]
    pub seed: u64,

    /// K-mer size
    #[arg(
        short = 'k',
        long,
        value_name = "INT",
        help_heading = "Sketching options"
    )]
    pub kmer: Option<u8>,
}

#[derive(Args, Debug)]
pub struct CompareArgs {}
