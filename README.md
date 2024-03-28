# darwin
[![Continuous Integration](https://github.com/Ebedthan/darwin/actions/workflows/ci.yml/badge.svg)](https://github.com/Ebedthan/darwin/actions/workflows/ci.yml)
<a href="https://github.com/Ebedthan/darwin/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue?style=flat">
</a>

darwin is a fast tool to build (rapid) neighbor-joining trees bases on mash distance. It takes as input the sequences (FASTA and FASTQ files are welcomed, compressed or not), compute the sketches and output a newick file of the tree.

The main advantages of darwin over others tools are:
- It uses the innovative approach of sketching algorithm [finch](https://github.com/onecodex/finch-rs) which is fast, have adaptive, count-based filtering (for FASTQs) and strandedness filtering. 
- Reliable and fast neighbor-joining tree estimation using [speedytree](https://docs.rs/speedytree/latest/speedytree/).

