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
[![DOI](https://zenodo.org/badge/243280413.svg)](https://zenodo.org/badge/latestdoi/243280413)

This project provides introductions to various bioinformatics tools with short guides, video demonstrations, and scripts that tie these tools together.
The documents in this project can be read locally in a plain-text editor, or viewed online at <https://rnnh.github.io/bioinfo-notebook/>.
If you are not familiar with using programs from the command line, begin with the page "[Introduction to the command line](docs/cl_intro.md)".
If you have any questions, or spot any mistakes, [please submit an issue on GitHub](https://github.com/rnnh/bioinfo-notebook/issues).

- [Pipeline examples](#pipeline-examples)
- [Contents](#contents)
- [Installation instructions](#installation-instructions)
- [Repository structure](#repository-structure)

## Pipeline examples

These bioinformatics pipelines can be carried out using scripts and tools described in this project.
Input files for some of these scripts can be specified in the command line; other scripts will need to be altered to fit the given input data.

### SNP analysis

- [FASTQ](docs/file_formats.md#fastq) reads from whole genome sequencing (WGS) can be assembled using [SPAdes](docs/SPAdes.md).
- Sequencing reads can be aligned to this assembled genome using [bowtie2](docs/bowtie2.md).
- The script [snp_calling.sh](docs/snp_calling.md) aligns sequencing reads to an assembled genome and detects single nucleotide polymorphisms (SNPs). This will produce a [Variant Call Format (VCF) file](docs/file_formats.md#vcf).
- The proteins in the assembled reference genome- the genome to which the reads are aligned- can be annotated using [genome_annotation_SwissProt_CDS.sh](docs/genome_annotation_SwissProt_CDS.md).
- The genome annotation [GFF](docs/file_formats.md#gff) file can be cross-referenced with the VCF file using [annotating_snps.R](docs/annotating_snps.md). This will produce an [annotated SNP format](docs/annotating_snps.md#annotated-snp-format) file.
- Annotated SNP format files can be cross-referenced using [annotated_snps_filter.R](docs/annotated_snps_filter.md). For two annotated SNP files, this script will produce a file with annotated SNPs unique to the first file, and a file with annotated SNPs unique to the second file.

### RNA-seq analysis

- [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md) can be used to download RNA-seq reads from NCBI's Sequence Read Archive (SRA) and align them to a reference genome. This script uses [fastq-dump](docs/fastq-dump.md) or [fasterq-dump](docs/fasterq-dump.md) to download the sequencing reads as [FASTQ](docs/file_formats.md#fastq), and [featureCounts](docs/featureCounts.md) to align them to a reference [FASTA nucleotide file.](docs/file_formats.md#fasta)
- Running [fastq-dump_to_featureCounts.sh](docs/fastq-dump_to_featureCounts.md) will produce feature count tables. These feature count tables can be combined using [combining_featCount_tables.py](docs/combining_featCount_tables.md).
- These combined feature count tables can be used for differential expression (DE) analysis. An example DE analysis script is included in this project: [DE_analysis_edgeR_script.R](docs/DE_analysis_edgeR_script.md). This script uses the [R programming language](https://cran.r-project.org/) with the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) library.

### Detecting orthologs between genomes

- [Augustus](docs/augustus.md) can be used to predict genes from [FASTA nucleotide files](docs/file_formats.md#fasta).
- Once the FASTA amino acid sequences have been [extracted from the Augustus annotations](docs/augustus.md#extracting-the-fasta-amino-acid-sequences-of-predicted-genes-from-an-augustus-annotation), you can search for orthologs using [OrthoFinder](docs/orthofinder.md).
- To find a specific gene of interest, search the amino acid sequences of the predicted genes using [BLAST](docs/blast.md).

## Contents

### [1. General guides](docs/part1.md)

- [Introduction to the command line](docs/cl_intro.md)
- [Windows Subsystem for Linux](docs/wsl.md)
- [Using Ubuntu through a Virtual Machine](docs/ubuntu_virtualbox.md)
- [File formats used in bioinformatics](docs/file_formats.md)

### [2. Program guides](docs/part2.md)

- [Augustus](docs/augustus.md)
- [Bcftools](docs/bcftools.md)
- [BLAST](docs/blast.md)
- [Bowtie](docs/bowtie.md)
- [Bowtie2](docs/bowtie2.md)
- [Conda](docs/conda.md)
- [Fasterq-dump](docs/fasterq-dump.md)
- [Fastq-dump](docs/fastq-dump.md)
- [FeatureCounts](docs/featureCounts.md)
- [Htseq-count](docs/htseq-count.md)
- [OrthoFinder](docs/orthofinder.md)
- [SAMtools](docs/samtools.md)
- [sgRNAcas9](docs/sgRNAcas9.md)
- [SPAdes](docs/SPAdes.md)

### [3. Scripts](docs/part3.md)

- [Annotated SNPs filter](docs/annotated_snps_filter.md)
- [Annotating SNPs](docs/annotating_snps.md)
- [Combining featCount tables.py](docs/combining_featCount_tables.md)
- [DE_analysis_edgeR_script.R](docs/DE_analysis_edgeR_script.md)
- [Fastq-dump to featureCounts](docs/fastq-dump_to_featureCounts.md)
- [Genome annotation script](docs/genome_annotation_SwissProt_CDS.md)
- [Linux setup script](docs/linux_setup.md)
- [SNP calling script](docs/snp_calling.md)
- [UniProt downloader](docs/UniProt_downloader.md)

## Installation instructions

After following these instructions, there will be a copy of the [bioinfo-notebook GitHub repo](https://www.github.com/rnnh/bioinfo-notebook/) on your system in the `~/bioinfo-notebook/` directory.
This means there will be a copy of all the documents and scripts in this project on your computer.
If you are using Linux and run the [Linux setup script](docs/linux_setup.sh), the `bioinfo-notebook` virtual environment- which includes the majority of the command line programs covered in this project- will also be installed using [conda](docs/conda.md).

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

**3.** If you are using Linux, run the [Linux setup script](docs/linux_setup.md) with this command after downloading the project:

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
│   ├── blastx_SwissProt_example_nucleotide_sequence.fasta.tsv
│   ├── blastx_SwissProt_S_cere.tsv
│   ├── design_table.csv
│   ├── example_genome_annotation.gtf
│   ├── example_nucleotide_sequence.fasta
│   └── featCounts_S_cere_20200331.csv
├── docs/
│   ├── annotated_snps_filter.md
│   ├── annotating_snps.md
│   ├── augustus.md
│   ├── blast.md
│   ├── bowtie2.md
│   ├── bowtie.md
│   ├── cl_intro.md
│   ├── cl_solutions.md
│   ├── combining_featCount_tables.md
│   ├── conda.md
│   ├── DE_analysis_edgeR_script.md
│   ├── DE_analysis_edgeR_script.pdf
│   ├── fasterq-dump.md
│   ├── fastq-dump.md
│   ├── fastq-dump_to_featureCounts.md
│   ├── featureCounts.md
│   ├── file_formats.md
│   ├── genome_annotation_SwissProt_CDS.md
│   ├── htseq-count.md
│   ├── linux_setup.md
│   ├── orthofinder.md
│   ├── part1.md    # Navigation page for website
│   ├── part2.md    # Navigation page for website
│   ├── part3.md    # Navigation page for website
│   ├── report_an_issue.md
│   ├── samtools.md
│   ├── sgRNAcas9.md
│   ├── snp_calling.md
│   ├── SPAdes.md
│   ├── ubuntu_virtualbox.md
│   ├── UniProt_downloader.md
│   └── wsl.md
├── envs/            # conda environment files
│   ├── augustus.yml            # environment for Augustus
│   ├── bioinfo-notebook.txt
│   ├── bioinfo-notebook.yml
│   ├── orthofinder.yml         # environment for OrthoFinder
│   └── sgRNAcas9.yml           # environment for sgRNAcas9
├── scripts/
│   ├── annotated_snps_filter.R
│   ├── annotating_snps.R
│   ├── combining_featCount_tables.py
│   ├── DE_analysis_edgeR_script.R
│   ├── fastq-dump_to_featureCounts.sh
│   ├── genome_annotation_SwissProt_CDS.sh
│   ├── linux_setup.sh
│   ├── snp_calling.sh
│   └── UniProt_downloader.sh
├── _config.yml     # Configures github.io project website
├── .gitignore
├── LICENSE
├── README.md
└── .travis.yml     # Configures Travis CI testing for GitHub repo
```
