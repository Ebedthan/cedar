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

    /// If some genomes can only be connected via divergent (statistically
    /// unreliable) pairs, search for a smaller dataset-wide k-mer size
    /// that resolves them, rather than refusing to build the tree.
    #[arg(long, help_heading = "Tree options")]
    pub include_div_pairs: bool,

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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn test_build_defaults() {
        let cli = Cli::try_parse_from(["cedar", "build", "genomes/"]).unwrap();
        match cli.command {
            Commands::Build(args) => {
                assert_eq!(args.indir, "genomes/");
                assert_eq!(args.sketch.size, 1000);
                assert_eq!(args.sketch.seed, 42);
                assert!(args.sketch.kmer.is_none());
                assert!(args.bootstrap.is_none());
                assert!(args.jacknife.is_none());
                assert_eq!(args.jacknife_prop, 0.5);
                assert!(!args.canonical);
                assert!(!args.include_div_pairs);
            }
            _ => panic!("expected Build"),
        }
    }

    #[test]
    fn test_dist_parses_with_output_and_sketch_overrides() {
        let cli = Cli::try_parse_from([
            "cedar", "dist", "genomes/", "-o", "out.tsv", "-s", "5000", "-k", "15",
        ])
        .unwrap();
        match cli.command {
            Commands::Dist(args) => {
                assert_eq!(args.indir, "genomes/");
                assert_eq!(args.output, Some("out.tsv".to_string()));
                assert_eq!(args.sketch.size, 5000);
                assert_eq!(args.sketch.kmer, Some(15));
            }
            _ => panic!("expected Dist"),
        }
    }

    #[test]
    fn test_bootstrap_and_jackknife_are_mutually_exclusive() {
        let result = Cli::try_parse_from(["cedar", "build", "genomes/", "-B", "100", "-J", "100"]);
        assert!(result.is_err(), "expected -B and -J to conflict");
    }

    #[test]
    fn test_bootstrap_alone_is_fine() {
        let cli = Cli::try_parse_from(["cedar", "build", "genomes/", "-B", "500"]).unwrap();
        match cli.command {
            Commands::Build(args) => assert_eq!(args.bootstrap, Some(500)),
            _ => panic!("expected Build"),
        }
    }

    #[test]
    fn test_dist_has_no_tree_or_bootstrap_flags() {
        // DistArgs should not accept build-only flags at all.
        let result = Cli::try_parse_from(["cedar", "dist", "genomes/", "-B", "100"]);
        assert!(
            result.is_err(),
            "dist should not accept -B (that's build-only)"
        );
    }

    #[test]
    fn test_missing_subcommand_errors_rather_than_defaulting() {
        let result = Cli::try_parse_from(["cedar"]);
        assert!(
            result.is_err(),
            "arg_required_else_help should reject no subcommand"
        );
    }
}
