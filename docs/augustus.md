---
layout: default
title: Augustus
parent: 2. Program guides
---

# Augustus

Augustus is a program that predicts genes in eukaryotic genomic sequences.
It can be run online, with a server for [smaller files](http://bioinf.uni-greifswald.de/augustus/submission.php) and one for [larger files](http://bioinf.uni-greifswald.de/webaugustus/), or locally.
The local version of Augustus can be installed through [conda](conda.md).
This project includes an example [augustus conda environment](../envs/augustus.yml).

## Predicting genes in a eukaryotic FASTA nucleic acid file using `augustus`

`augustus` can be used to predict genes as follows:

```bash
$ augustus --species=species_name input_file.fna > output_file.gff
```

In this command...

1. `--species` is used to specify the target species for gene predictions (`species_name`).
2. `input_file.fna` is the input FASTA nucleic acid file ([.fna](file_formats.md#fasta)).
3. `output_file.gff` is the general feature format ([GFF](file_formats.md#generic-feature-formats)) genome annotation output file.
Lines beginning with `#` are Augustus comments: these lines do not follow the GFF structure.

The following command gives the list of valid species names for use with Augustus:

```bash
$ augustus --species=help
```

## Extracting the FASTA amino acid sequences of predicted genes from an Augustus annotation

The genome annotation file produced by `augustus` (`output_file.gff`) contains the amino acid sequences of predicted genes in comment lines.
These amino acid sequences can be extracted to a FASTA file with the following command:

```bash
$ getAnnoFasta.pl output_file.gff
```

The amino acid sequences will be written to `output_file.aa`.
This is a FASTA amino acid ([.faa](file_formats.md#fasta)).
The extension of this file can be changed from ".aa" to ".faa" with the following command:

```bash
$ mv output_file.aa output_file.faa
```

## Removing comments from Augustus annotations

Genome annotations produced by Augustus follow the [Generic Feature Format](file_formats.md#generic-feature-formats), with the addition of comment lines for amino acid sequences.
These are the same FASTA amino acid sequences that are extracted using `getAnnoFasta.pl`.
These lines begin with the character `#`, and removing them results a standard GFF file.

Here is one method for removing these amino acid lines, using `grep -v` to select lines which do not contain the `#` character:

```bash
$ grep -v "#" augustus_annotation.gff > clean_augustus_annotation.gff
```

## Demonstration

In this video, `augustus` is used to predict genes in `example_nucleotide_sequence.fasta`.
This results in a genome annotation file: `augustus_example.gff`.
The script `getAnnoFasta.pl` is used to extract the amino acid sequences in this genome annotation file to a new FASTA amino acid file: `augustus_example.aa`.
The `mv` command is used to change the extension of this file from ".aa" to ".faa".

[![asciicast](https://asciinema.org/a/346541.svg)](https://asciinema.org/a/346541?autoplay=1)

## See also

- [conda](conda.md)
- [augustus conda environment](../envs/augustus.yml)
- [File formats used in bioinformatics](file_formats.md)

## References

- [The Augustus website](http://bioinf.uni-greifswald.de/augustus/)
- [GNU grep](https://www.gnu.org/software/grep/manual/grep.html)
