# darwin
[![Continuous Integration](https://github.com/Ebedthan/darwin/actions/workflows/ci.yml/badge.svg)](https://github.com/Ebedthan/darwin/actions/workflows/ci.yml)
<a href="https://github.com/Ebedthan/darwin/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue?style=flat">
</a>

darwin is a fast tool to build (rapid) neighbor-joining trees bases on mash distance. It takes as input the sequences (FASTA and FASTQ files are welcomed, compressed or not), compute the sketches and output a newick file of the tree.

The main advantages of darwin over others tools are:
- It uses the innovative approach of sketching algorithm [finch](https://github.com/onecodex/finch-rs) which is fast, have adaptive, count-based filtering (for FASTQs) and strandedness filtering. 
- Reliable and fast neighbor-joining tree estimation using [speedytree](https://docs.rs/speedytree/latest/speedytree/).

darwin outputs the tree in newick format.

# Quick start guide

```
# Compute rapid neighbor-joining tree of all files in a directory
darwin dir/*

# Compute rapid NJ tree using specific files
darwin file1.fa.gz file2.fq.xz file3.fna.bz2

# Compute canonical neighbor-joining tree
darwin -c dir/*
```
Full help is available from `darwin --help`;

# Installation

```
git clone https://github.com/Ebedthan/darwin.git
cd darwin

# if default Rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
darwin -h
```
### Minimum supported Rust version
`darwin` minimum [Rust](https://www.rust-lang.org/) version is 1.74.1.

### Semver
`darwin` is following [Semantic Versioning 2.0](https://semver.org/).

### Licence
`darwin` is distributed under the terms of the MIT license.
See [LICENSE-MIT](https://github.com/Ebedthan/xgt/blob/main/LICENSE-MIT) for details.