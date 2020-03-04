# README

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project aims to give brief introductions for various bioinformatics tools with quick start guides. If you know how to use the command line, hopefully these guides will get you started with each program in a matter of minutes. These guides also include video demonstrations, created using [Asciinema](https://asciinema.org/~rnnh). This site will also feature scripts that tie together the various programs covered.

If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Quick Start Guides

1. [bowtie2](docs/bowtie2.md) *For aligning reads to sequences*
2. [conda](docs/conda.md) *For downloading and sharing sets of programs*
3. [fastq-dump](docs/fastq-dump.md) *For downloading reads from NCBI's SRA*
4. [featureCounts](docs/featureCounts.md) *For assigning alignments to genome annotation features*
5. [htseq-count](docs/htseq-count.md) *For aligning reads to genomic features*
6. [samtools](docs/samtools.md) *For manipulating alignment files*

## Scripts

1. [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md). This script downloads FASTQ reads, aligns them to a reference genome, and counts how many reads align to each gene in that reference genome.
