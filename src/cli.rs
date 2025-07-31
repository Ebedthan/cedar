// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "darwin",
    about = "Compute (rapid) neighbor joining tree from sequences",
    author,
    version,
    arg_required_else_help = true
)]
pub struct Cli {
    /// Fasta file(s) to build trees [suports .gz, .xz, .bz2]
    #[arg(required = true)]
    pub input: Vec<String>,

    /// Output tree (Newick format) to FILE
    #[arg(short, value_name = "FILE")]
    pub output: Option<String>,

    /// Keep sketches and distance files
    #[arg(short = 'K')]
    pub keep: bool,

    /// Number of threads to use
    #[arg(short, default_value_t = 1, value_name = "INT")]
    pub threads: usize,

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
    pub kmer: u8,

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
