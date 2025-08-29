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
    /// Build NJ tree using mash distance
    Build(BuildArgs),

    /// Compare trees
    Compare(CompareArgs),
}
#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Fasta file(s) to build trees [suports .gz, .xz, .bz2]
    #[arg(required = true)]
    pub input: Vec<String>,

    /// Output tree (Newick format) to FILE
    #[arg(short, value_name = "FILE")]
    pub output: Option<String>,

    /// Keep sketches and distance files
    #[arg(short = 'K')]
    pub keep: bool,

    /// Temporary directory name
    #[arg(long, value_name = "DIR", default_value_t = String::from("cedar_tmp"))]
    pub tempdir: String,

    /// Boostrap replicates
    #[arg(short, value_name = "INT")]
    pub bootstrap: Option<usize>,

    /// Sketch size
    #[arg(
        short = 's',
        long,
        default_value_t = 1000,
        value_name = "INT",
        help_heading = "Sketching options"
    )]
    pub size: usize,

    /// Seed for the hash function
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

    /// Amount of extra scketching before filtering
    #[arg(
        short = 'x',
        long,
        default_value_t = 200,
        value_name = "INT",
        help_heading = "Sketching options"
    )]
    pub oversketch: usize,

    /// Compute canonical NJ tree
    #[arg(short = 'c', help_heading = "Tree options")]
    pub canonical: bool,
}

#[derive(Args, Debug)]
pub struct CompareArgs {}
