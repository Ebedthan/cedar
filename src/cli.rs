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
    #[arg(short, default_value_t = 1, value_name = "INT")]
    pub threads: usize,

    #[arg(short = 'v', long, action = clap::ArgAction::SetTrue)]
    pub verbose: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Build a Mash-distance neighbor-joining tree
    Build(BuildArgs),
    // `Compare` intentionally omitted until it's implemented — see roadmap.
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
