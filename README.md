# README

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project aims to give brief introductions for various bioinformatics tools with quick start guides, and provide scripts that tie these tools together. If you know how to use the command line, hopefully these guides will get you started with each program in a matter of minutes.

For reproducibility, a [conda](docs/conda.md) environment is used throughout this project. Video demonstrations created using [Asciinema](https://asciinema.org/~rnnh) are also provided.

If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Quick Start Guides

1. [bowtie2](docs/bowtie2.md) *aligns reads to sequences*
2. [conda](docs/conda.md) *facilitates downloading and sharing sets of programs*
3. [fastq-dump](docs/fastq-dump.md) *downloads reads from NCBI's SRA*
4. [featureCounts](docs/featureCounts.md) *assigns alignments to genome annotation features*
5. [htseq-count](docs/htseq-count.md) *aligns reads to genomic features*
6. [samtools](docs/samtools.md) *manipulates alignment files and simluates sequence reads*

## Scripts

1. [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md). This script downloads FASTQ reads, aligns them to a reference genome, and counts how many reads align to each gene in that reference genome.
