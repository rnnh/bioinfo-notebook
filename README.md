---
layout: default
title: Home
nav_order: 1
description: "Quick start guides for bioinformatics programs, with video demonstrations and scripts."
permalink: /
---


# [Bioinformatics Notebook](https://github.com/rnnh/bioinfo-notebook.git)

by [Ronan Harrington](https://github.com/rnnh)

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project provides introductions to various bioinformatics tools with short guides, video demonstrations, and scripts that tie these tools together.
The documents in this project can be read locally in a plain-text editor, or viewed online at <https://rnnh.github.io/bioinfo-notebook/>.

If you are not familiar with using programs from the command line, begin with the page "[Introduction to the command line](docs/cl_intro.md)".

If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Contents

- [Installation instructions](#installation-instructions)
	- [Video demonstration of installation](#Video-demonstration-of-installation)
- [Repository structure](#repository-structure)

### [1. General guides](docs/part1.md)

- [Introduction to the command line](docs/cl_intro.md)
- [Windows Subsystem for Linux](docs/wsl.md)
- [Using Ubuntu through a Virtual Machine](docs/ubuntu_virtualbox.md)
- [File formats used in bioinformatics](docs/file_formats.md)

### [2. Program guides](docs/part2.md)

- [Augustus](docs/augustus.md)
- [BLAST](docs/blast.md)
- [bowtie](docs/bowtie.md)
- [bowtie2](docs/bowtie2.md)
- [conda](docs/conda.md)
- [fasterq-dump](docs/fasterq-dump.md)
- [fastq-dump](docs/fastq-dump.md)
- [featureCounts](docs/featureCounts.md)
- [htseq-count](docs/htseq-count.md)
- [samtools](docs/samtools.md)

### [3. Scripts](docs/part3.md)

- [linux_setup.sh](docs/linux_setup.md)
- [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md)
- [combining_featCount_tables.py](docs/combining_featCount_tables.md)

## Installation instructions

After following these instructions, there will be a copy of the [bioinfo-notebook GitHub repo](https://www.github.com/rnnh/bioinfo-notebook/) on your system in the `~/bioinfo-notebook/` directory.
This means there will be a copy of all the documents and scripts in this project on your computer.
The `bioinfo-notebook` virtual environment, which includes all of the command line programs covered in this project, will also be installed on your Ubuntu system using [conda](docs/conda.md).

**1.** This project is written to be used through a UNIX (Linux or Mac with macOS Mojave or later) operating system.
 If you are using a Windows operating system, begin with these pages on setting up Ubuntu (a Linux operating system):
 
- [Windows Subsystem for Linux](docs/wsl.md)
- [Using Ubuntu through a Virtual Machine](docs/ubuntu_virtualbox.md)

Once you have an Ubuntu system set up, run the following command to update the lists of available software:

```bash
$ sudo apt-get update # Updates lists of software that can be installed
```

**2.** Run the following command in your home directory (`~`) to download this project:

```bash
$ git clone https://github.com/rnnh/bioinfo-notebook.git
```

**3.** After downloading this project, run the [Linux setup script](docs/linux_setup.md) with this command:

```bash
$ bash ~/bioinfo-notebook/scripts/linux_setup.sh
```

### Video demonstration of installation

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)

## Repository structure

```
bioinfo-notebook/
├── assets/
│   └── bioinfo-notebook_logo.svg
├── data/
│   ├── example_genome_annotation.gtf
│   └── example_nucleotide_sequence.fasta
├── docs/
│   ├── Augustus.md
│   ├── blast.md
│   ├── bowtie.md
│   ├── bowtie2.md
│   ├── cl_intro.md
│   ├── cl_solutions.md
│   ├── combining_featCount_tables.md
│   ├── conda.md
│   ├── fasterq-dump.md
│   ├── fastq-dump.md
│   ├── fastq-dump_to_featureCounts.md
│   ├── featureCounts.md
│   ├── file_formats.md
│   ├── htseq-count.md
│   ├── linux_setup.md
│   ├── part1.md    # Navigation page for website
│   ├── part2.md    # Navigation page for website
│   ├── part3.md    # Navigation page for website
│   ├── samtools.md
│   ├── to_do.md
│   ├── ubuntu_virtualbox.md
│   └── wsl.md
├── envs/            # conda environment files
│   ├── augustus.txt # environments for augustus
│   ├── augustus.yml
│   ├── bioinfo-notebook.txt
│   └── bioinfo-notebook.yml
├── scripts/
│   ├── combining_featCount_tables.py
│   ├── fastq-dump_to_featureCounts.sh
│   └── linux_setup.sh
├── _config.yml     # Configures github.io project website
├── .gitignore
├── LICENSE
├── README.md
└── .travis.yml     # Configures Travis CI testing
```
