// Copyright 2024-2025 Anicet Ebou.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed except according
// to those terms.

use clap::{crate_authors, crate_version, Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    Command::new("darwin")
        .about("darwin - compute (rapid) neighbor joining tree from sequences")
        .author(crate_authors!())
        .version(crate_version!())
        .arg_required_else_help(true)
        .args([
            Arg::new("INPUT")
                .help("The file(s) to build trees [supports .gz, .xz, .bz2]")
                .num_args(1..)
                .required(true),
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output tree (Newick format) to FILE"),
            Arg::new("keep")
                .short('K')
                .long("keep")
                .help("Keep sketches and distance files")
                .action(ArgAction::SetTrue),
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("INT")
                .default_value("1")
                .help("Number of threads to use"),
        ])
        .args(build_sketching_args())
        .args(build_tree_args())
}

/// Helper function to create sketching-related arguments
fn build_sketching_args() -> [Arg; 4] {
    [
        Arg::new("size")
            .short('s')
            .long("size")
            .value_name("INT")
            .default_value("1000")
            .help("Sketch size")
            .help_heading("Sketching options"),
        Arg::new("seed")
            .short('S')
            .long("seed")
            .value_name("INT")
            .default_value("42")
            .help("Seed for the hash function")
            .help_heading("Sketching options"),
        Arg::new("kmer")
            .short('k')
            .long("kmer")
            .value_name("INT")
            .default_value("21")
            .help("K-mer size")
            .help_heading("Sketching options"),
        Arg::new("oversketch")
            .short('x')
            .long("oversketch")
            .value_name("INT")
            .default_value("200")
            .help("Amount of extra sketching before filtering")
            .help_heading("Sketching options"),
    ]
}

/// Helper function to create tree-related arguments
fn build_tree_args() -> [Arg; 1] {
    [Arg::new("canonical")
        .short('c')
        .long("canonical")
        .help("Compute canonical NJ tree")
        .action(ArgAction::SetTrue)
        .help_heading("Tree options")]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_cmd() {
        build_cli().debug_assert();
    }
}
