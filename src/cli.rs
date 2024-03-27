use clap::{crate_authors, crate_version, Arg, ArgAction, Command};

pub fn build_cli() -> Command {
    Command::new("darwin")
        .about("darwin - compute (rapid) neighbor joining tree from sequences")
        .arg_required_else_help(true)
        .author(crate_authors!())
        .version(crate_version!())
        .arg(
            Arg::new("INPUT")
                .help("the file(s) to build trees [can be .gz, .xz, .bz2]")
                .num_args(1..)
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .value_name("FILE")
                .help("output tree (newick format) to FILE"),
        )
        .arg(
            Arg::new("keep")
                .short('K')
                .help("keep sketches and distances files")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .value_name("INT")
                .help("number of threads")
                .default_value("1"),
        )
        .arg(
            Arg::new("size")
                .short('s')
                .value_name("INT")
                .help("sketch size")
                .default_value("1000")
                .help_heading("Sketching options"),
        )
        .arg(
            Arg::new("seed")
                .short('S')
                .value_name("INT")
                .help("seed to provide to the hash function")
                .default_value("42")
                .help_heading("Sketching options"),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .value_name("INT")
                .default_value("21")
                .help("k-mer size")
                .help_heading("Sketching options"),
        )
        .arg(
            Arg::new("oversketch")
                .short('x')
                .value_name("INT")
                .help("amount of extra sketching to do before filtering")
                .default_value("200")
                .help_heading("Sketching options"),
        )
        .arg(
            Arg::new("canonical")
                .short('c')
                .action(ArgAction::SetTrue)
                .help("compute canonical NJ tree")
                .help_heading("Tree options"),
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_cmd() {
        build_cli().debug_assert();
    }
}
