---
layout: default
title: BLAST
parent: 2. Program guides
---

# BLAST

The Basic Local Alignment Search Tool (BLAST) is an algorithm and program for comparing primary biological sequence information, such as the amino-acid sequences of proteins or the nucleotides of DNA and/or RNA sequences.
BLAST is one of the most widely used tools in bioinformatics; it can be applied to different problems or projects in a myriad ways.

BLAST can be used online, or through the command line.
Most biologists are familiar with NCBI's web application for BLAST (<https://blast.ncbi.nlm.nih.gov/Blast.cgi>).
If you use this web application regularly, the command line BLAST program is worth your consideration.
The command line version of BLAST has several advantages over the web version:

1. BLAST on the command line can be used to run *local searches*, i.e. searches which use files that are on your computer, instead of files that are on an NCBI database.
2. BLAST searches on the command line can be made more specific by adding additional arguments.
3. BLAST searches carried out on the command line can be automated, and incorporated into larger scripts.
4. The command line BLAST program can output search results in various structured text formats.

## Creating a BLAST database using `makeblastdb`

To search against a set of nucleotide or amino acid sequences using BLAST, a database must be created.
This can be done using the `makeblastdb` command.

```bash
$ makeblastdb -dbtype prot/nucl -in input_file -out database_name
```

In this command...

1. `-dbtype` specifies the type of sequences used to create the database. For amino acid (protein) sequences, `prot` is used ("`-dbtype prot`"). For nucleic acid sequences, `nucl` is used ("`-dbtype nucl`").
2. `-in` is used to specify the input file. The database created can be used to search against the sequences in this file.
3. `-out` is used to name the database that will be created from the input file.

## Video demonstration

In this demonstration, `makeblastdb` is used to create a BLAST database from the file `S_cere_genomes.fna`.
This FASTA nucleic acids (`.fna`) file was created by concatenating the following *Saccharomyces cerevisiae* genome assemblies, which were downloaded from NCBI: [GCA_003086655.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_003086655.1), [GCA_003709285.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_003709285.1) and [GCA_004328465.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_004328465.1).

The program `blastn` is then used to query `23S_rRNA_gene.fna` against this database.
This file is a copy of the [*Scutellospora fulgida* isolate NC303A 25S ribosomal RNA gene](https://www.ncbi.nlm.nih.gov/nuccore/FJ461870.1?report=fasta) from NCBI.

The program `tblastn` is also used to query `YPK1.faa` against this database multiple times.
This FASTA amino acid (`.faa`) file is a copy of the [serine/threonine-protein kinase YPK1](https://www.uniprot.org/uniprot/P12688) from UniProt.
This search is carried out multiple times with additional parameters: the flag `-evalue`is used to set an E-value threshold, and the flag `-max_hsps` is used to set a maximum number of High-scoring Segment Pairs (HSPs).

The results from these BLAST searches are written to tab-separated values (`.tsv`) files.
This output format is specified with the flag `-outfmt 6`.

[![asciicast](https://asciinema.org/a/327279.svg)](https://asciinema.org/a/327279?autoplay=1)
