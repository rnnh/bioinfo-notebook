# [Bioinformatics Notebook](https://rnnh.github.io/bioinfo-notebook/)

by [Ronan Harrington](https://github.com/rnnh)

[![Build Status](https://travis-ci.com/rnnh/bioinfo-notebook.svg?branch=master)](https://travis-ci.com/rnnh/bioinfo-notebook)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/rnnh/bioinfo-notebook)
![GitHub repo size](https://img.shields.io/github/repo-size/rnnh/bioinfo-notebook)
![Website](https://img.shields.io/website?url=https%3A%2F%2Frnnh.github.io%2Fbioinfo-notebook)

This project provides introductions to various bioinformatics tools with short guides, video demonstrations, and scripts that tie these tools together.
The documents in this project can be read locally in a plain-text editor, or viewed online at <https://rnnh.github.io/bioinfo-notebook/>.

If you have any suggestions, or questions, or spot any mistakes, [please let me know](https://github.com/rnnh/bioinfo-notebook/issues).

## Contents

- [Installation instructions](#installation-instructions)
	- [Video demonstration of installation](#Video-demonstration-of-installation)

### Part I: General guides

- [Windows Subsystem for Linux (Windows only)](docs/wsl.md)
- [Using Ubuntu through a Virtual Machine (Mac or Windows)](docs/ubuntu_virtualbox.md)
- [File formats used in bioinformatics](docs/file_formats.md)

### Part II: Program guides

- [bowtie2](docs/bowtie2.md)
- [conda](docs/conda.md)
- [fastq-dump](docs/fastq-dump.md)
- [featureCounts](docs/featureCounts.md)
- [htseq-count](docs/htseq-count.md)
- [samtools](docs/samtools.md)

### Part III: Scripts

- [linux_setup.sh](docs/linux_setup.md)
- [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md)
- [combining_featCount_tables.py](docs/combining_featCount_tables.md)

## Installation instructions

1. This project is written to be used through the Ubuntu operating system.
 If you are using a Windows/Mac operating system, begin with these pages on setting up Ubuntu:
	- [Windows Subsystem for Linux (Windows only)](docs/wsl.md)
	- [Using Ubuntu through a Virtual Machine (Mac or Windows)](docs/ubuntu_virtualbox.md)

2. Once you have an Ubuntu system set up, run the following command to update the lists of available software:

```bash
$ sudo apt-get update # Updates lists of software that can be installed
```

2. Run the following command in your home directory (`~`) to download this project:

```bash
$ git clone https://github.com/rnnh/bioinfo-notebook.git
```

3. After downloading this project, run the [Linux setup script](docs/linux_setup.md) with this command:

```bash
$ bash ~/bioinfo-notebook/scripts/linux_setup.sh
```

After following these instructions, there will be a copy of the [bioinfo-notebook GitHub repo](https://www.github.com/rnnh/bioinfo-notebook/) on your Ubuntu system in the `~/bioinfo-notebook/` directory.
This means there will be a copy of all the documents and scripts in this project on your computer.
The `bioinfo-notebook` virtual environment, which includes all of the command line programs covered in this project, will also be installed on your Ubuntu system using [conda](docs/conda.md).

### Video demonstration of installation

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)
