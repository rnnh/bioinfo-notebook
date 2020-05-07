---
layout: default
title: BLAST
parent: 2. Program guides
---

# BLAST

The Basic Local Alignment Search Tool (BLAST) is an algorithm and program for comparing primary biological sequence information, such as the amino-acid sequences of proteins or the nucleotides of DNA and/or RNA sequences.
BLAST is one of the most widely used tools in bioinformatics; it can be applied to different problems or projects in a myriad ways.

BLAST can be used online, or through the command line.
Most biologists are familiar with [NCBI's web application for BLAST](<https://blast.ncbi.nlm.nih.gov/Blast.cgi>).
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


## Searching against a BLAST nucleotide database using `blastn`

The program `blastn` is used for searching nucleotide databases with a nucleotide query.

```bash
$ blastn -query query_file.fna -db nucl_database_name -out results_file.tsv -outfmt 6 -evalue x -max_hsps y
```

In this command...

1. `-query` is used to select the FASTA nucleic acids file you want to search against the BLAST database (the `query_file.fna`).
2. `-db` is used to select the BLAST nucleotide database you want to search against (`nucl_database_name`).
3. `-out` is used to direct the results to an output file (`results_file.tsv`).
4. `-outfmt` is used to specify how this results file should be formatted. In this case, as `-outfmt` is `6`, the results will be written to a file as tab-separated values: this is why `results_file.tsv` has a `.tsv` extension.
5. `-evalue` is used to set an E-value threshold (`x`). Results which have an E-value greater than this threshold will not be written to the results file.
6. `-max_hsps` is used to set a High-scoring Segment Pairs (HSPs) threshold (`y`). When given, no more than `y` HSPs (alignments) for each query-subject pair will be written to the results file.

The last two arguments given in this command- `-evalue` and `-max_hsps`- are optional, but they are useful as they allow the results to be filtered before being written to the file.
Using these arguments will result in more specific results, and will reduce the need to manually filter results later.

## BLAST `-outfmt 6` results

These BLAST results are taken from the [video demonstration](#video_demonstration) and are in BLAST output format 6.

```
gi|242120357|gb|FJ461870.1|	NC_001144.5	93.252	163	11	0	196	358	454921	454759	7.57e-63	241
gi|242120357|gb|FJ461870.1|	NC_001144.5	93.252	163	11	0	196	358	464058	463896	7.57e-63	241
gi|242120357|gb|FJ461870.1|	CP036478.1	93.252	163	11	0	196	358	454829	454667	7.57e-63	241
gi|242120357|gb|FJ461870.1|	CP036478.1	93.252	163	11	0	196	358	463966	463804	7.57e-63	241
gi|242120357|gb|FJ461870.1|	CP024006.1	93.252	163	11	0	196	358	453978	453816	7.57e-63	241
```

These results are tab-separated values, meaning each column in the results is separated by a `Tab` character.
These columns always appear in the same order:

```
query_id	subject_id	per_identity	aln_length	mismatches	gap_openings	q_start	q_end	s_start	s_end	e-value	bit_score
```

In this format...

1. `query_id` is the FASTA header of the sequence being searched against the database (the query sequence).
2. `subject_id` is the FASTA header of the sequence in the database that the query sequence has been aligned to (the subject sequence).
3. `per_identity` is the percentage identity- the extent to which the query and subject sequences have the same residues at the same positions.
4. `aln_length` is the alignment length.
5. `mismatches` is the number of mismatches.
6. `gap_openings` is the number of gap openings in the alignment.
7. `q_start` is the start of the alignment in the query sequence.
8. `q_end` is the end of the alignment in the query sequence.
9. `s_start` is the start of the alignment in the subject sequence.
10. `s_end` is the end of the alignment in the subject sequence.
11. `e_value` is the expect value (E-value) for the alignment.
12. `bit_score` is the bit-score of the alignment.

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
