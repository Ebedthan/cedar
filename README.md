# cedar: from genomes to trees, simplified.

[![CI](https://github.com/Ebedthan/cedar/actions/workflows/ci.yml/badge.svg)](https://github.com/Ebedthan/cedar/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/Ebedthan/cedar/graph/badge.svg?token=S3OLFFHF4X)](https://codecov.io/gh/Ebedthan/cedar)
[![License](https://img.shields.io/badge/license-MIT-blue?style=flat)](https://github.com/Ebedthan/cedar/blob/main/LICENSE)


## 🗺️ Overview

`cedar` is a fast, portable, and reproducible toolkit for phylogenomics, written in pure Rust.
It streamlines the reconstruction and comparison of species trees from whole genomes or sets of orthologous groups.
With sensible defaults, parallel execution, and a clean CLI, Cedar is designed to be as easy to use for biologists as it is powerful for computational phylogeneticists.

## 🌲 Cedar Roadmap
### ✅ Core v1.0

- **Input handling**
  - Accept genome FASTA files (.fa/.fna/.faa, compressed or uncompressed).
  - Accept orthologous group FASTAs (OMA, OrthoFinder, BUSCO outputs).
  - Transparent handling of compressed inputs (.gz, .bz2, .xz).
- **Supermatrix workflow**
  - Multiple sequence alignment (MAFFT wrapper).
  - Concatenate per-locus alignments into a supermatrix.
  - Generate partition files for downstream ML tree inference.
  - Handle missing data (gap padding, min. occupancy thresholds).
- **Tree building**
  - Distance-based tree inference using Mash.
  - Wrapper for IQ-TREE and FastTree.
  - Support user-defined model options.
- **Tree comparison**
  - Robinson-Foulds distance.
  - Leaf reconciliation (ignore taxa missing in one tree).
- **CLI design**
  - Subcommands:
    - cedar build → build species trees.
    - cedar compare → compare two trees.
  - Smart defaults (auto-detect dataset size for MAFFT mode, etc.).
  - Clear progress bars and logging.
- **Reproducibility**
  - Output manifest (YAML/JSON) with command-line options, versions, and seeds.
  - Deterministic runs by default.

### Minimum supported Rust version
`cedar` minimum [Rust](https://www.rust-lang.org/) version is 1.74.1.

### Semver
`cedar` is following [Semantic Versioning 2.0](https://semver.org/).

### Licence
`cedar` is distributed under the terms of the MIT license.
See [LICENSE-MIT](https://github.com/Ebedthan/cedar/blob/main/LICENSE-MIT) for details.
