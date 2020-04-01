# Bioinformatics Notebook

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project aims to give brief introductions for various bioinformatics tools with quick start guides, and provide scripts that tie these tools together.
For reproducibility, [conda](docs/conda.md) environments are used throughout this project.
Video demonstrations created using [Asciinema](https://asciinema.org/~rnnh) are also provided.
If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Installation

1. This project is written to be used through the Ubuntu operating system.
 If you are using a Windows/Mac operating system, begin with these pages on setting up Ubuntu:
	- [Windows Subsystem for Linux (Windows only)](docs/wsl.md)
	- [Using Ubuntu through a Virtual Machine (Mac or Windows)](docs/ubuntu_virtualbox.md)

2. Once you have an Ubuntu system set up, run the following command in your home directory (`~`) to download this project:

```bash
$ git clone https://github.com/rnnh/bioinfo-notebook.git
```

3. After downloading this project, run the [Linux setup script](docs/linux_setup.md) with these commands:

```bash
$ sudo apt-get update # Updates lists of software that can be installed
$ bash ~/bioinfo-notebook/scripts/linux_setup.sh
```

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
