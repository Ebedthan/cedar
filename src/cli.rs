// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
    name = "cedar",
    about = "Uncertainty-aware distance-based phylogenomics from Mash sketches",
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
    /// Build a Mash-distance neighbor-joining consensus tree
    Build(BuildArgs),

    /// Compute pairwise Mash distances with uncertainty estimates (no tree)
    Dist(DistArgs),
}

/// Sketching options shared between `build` and `dist`: either set sketch
/// size / k-mer size manually, or let cedar choose them to hit a target
/// precision on each pairwise distance estimate.
#[derive(Args, Debug)]
pub struct SketchArgs {
    /// Sketch size
    #[arg(
        short = 's',
        long,
        default_value_t = 1000,
        value_name = "INT",
        conflicts_with = "target_precision",
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
        conflicts_with = "target_precision",
        help_heading = "Sketching options"
    )]
    pub kmer: Option<u8>,

    /// Automatically choose sketch size (and k-mer size) so that each
    /// pairwise Mash distance's relative uncertainty stays under this
    /// fraction (e.g. 0.1 for ~10%), instead of a fixed sketch size.
    #[arg(
        long,
        value_name = "FLOAT",
        conflicts_with_all = ["size", "kmer"],
        help_heading = "Sketching options"
    )]
    pub target_precision: Option<f64>,
}

#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Input directory containing FASTA files
    #[arg(value_name = "DIR")]
    pub indir: String,

    /// Output tree (Newick format) to FILE
    #[arg(short, value_name = "FILE")]
    pub output: Option<String>,

    /// Bootstrap replicates
    #[arg(
        short = 'B',
        long = "boot",
        value_name = "INT",
        conflicts_with = "jacknife",
        help_heading = "Tree options"
    )]
    pub bootstrap: Option<usize>,

    /// Jackknife replicates
    #[arg(
        short = 'J',
        long = "jack",
        value_name = "INT",
        conflicts_with = "bootstrap",
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

    #[command(flatten)]
    pub sketch: SketchArgs,
}

#[derive(Args, Debug)]
pub struct DistArgs {
    /// Input directory containing FASTA files
    #[arg(value_name = "DIR")]
    pub indir: String,

    /// Output distance table (TSV) to FILE; prints to stdout if omitted
    #[arg(short, value_name = "FILE")]
    pub output: Option<String>,

    #[command(flatten)]
    pub sketch: SketchArgs,
}
