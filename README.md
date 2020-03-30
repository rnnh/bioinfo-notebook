# README

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project aims to give brief introductions for various bioinformatics tools with quick start guides, and provide scripts that tie these tools together.

If you would like to try the scripts or the commands in the guides on this website, a Linux system such as Ubuntu is recommended.
If you are using a Windows or Mac computer, you can use Ubuntu through [VirtualBox](docs/ubuntu_virtualbox.md) or [Windows Subsystem for Linux](docs/wsl.md).

For reproducibility, a [conda](docs/conda.md) environment is used throughout this project. Video demonstrations created using [Asciinema](https://asciinema.org/~rnnh) are also provided.

If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Quick Start Guides

- [File formats](docs/file_formats.md) *brief introduction to bioinformatics file formats*

### Using Ubuntu through a Windows or Mac computer

- [Using Ubuntu through a Virtual Machine](docs/ubuntu_virtualbox.md)
- [Windows Subsystem for Linux](docs/wsl.md)

### Command line programs

- [bowtie2](docs/bowtie2.md) *aligns reads to sequences*
- [conda](docs/conda.md) *facilitates downloading and sharing sets of programs*
- [fastq-dump](docs/fastq-dump.md) *downloads reads from NCBI's SRA*
- [featureCounts](docs/featureCounts.md) *assigns alignments to genome annotation features*
- [htseq-count](docs/htseq-count.md) *aligns reads to genomic features*
- [SAMtools](docs/samtools.md) *manipulates alignment files and simluates sequence reads*

## Scripts

- [linux_setup.sh](docs/linux_setup.md). This script downloads and installs Miniconda3, and installs the bioinfo-notebook virtual environment using conda. Video demonstration:

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)

- [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md). This script downloads FASTQ reads, aligns them to a reference genome, and counts how many reads align to each gene in that reference genome. Video demonstration:

[![asciicast](https://asciinema.org/a/308745.svg)](https://asciinema.org/a/308745?autoplay=1)

- [combining_featCount_tables.py](docs/combining_featCount_tables.md) This is a Python script that creates a single CSV feature count table from the featureCounts output tables in the target directory. Video demonstration:

[![asciicast](https://asciinema.org/a/311771.svg)](https://asciinema.org/a/311771?autoplay=1)
