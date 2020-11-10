---
layout: default
title: sgRNAcas9
parent: 2. Program guides
---

# sgRNAcas9

sgRNAcas9 is [a package for designing CRISPR sgRNA and evaluating potential off-target cleavage sites](https://doi.org/10.1371/journal.pone.0100448).

## Running sgRNAcas9

1. Install the [conda](conda.md) virutal environment for [sgRNAcas9](../envs/sgRNAcas9.yml).
2. Download [the GUI version of sgRNAcas9 from SourceForge](https://sourceforge.net/projects/sgrnacas9/).
3. Activate the sgRNAcas9 virtual environment.
4. In the directory for sgRNAcas9, run the following command to launch the sgRNAcas9 graphical user interface (GUI):

```bash
(sgRNAcas9) ~/sgRNAcas9_V3.0_GUI$ java -jar sgRNAcas9.jar
```

## Using sgRNAcas9

In the sgRNAcas9 GUI...

- Select the [FASTA nucleic acid](file_formats.md#fasta) file of the target sequences in the "Target sequences(FASTA):" dialog box.
- Select the [FASTA nucleic acid](file_formats.md#fasta) file of the genome you want to design the guide RNAs for in the "Genome sequence(FASTA):" dialog box.
- Click "RUN" to run the program

sgRNAcas9 will create a `report` directory in the current working directory.
This directory contains its results.
The most important file in this directory is `sgRNAcas9_report.xls`.
This Excel files contains reported guide RNA sequences for CRISPR with quality score, and counts of potential off-target sites.

## References

- [sgRNAcas9 paper](https://sourceforge.net/projects/sgrnacas9/)
- [sgRNAcas9 website](http://biootools.com/software.html)
- [sgRNAcas9 on SourceForge](https://sourceforge.net/projects/sgrnacas9/)
