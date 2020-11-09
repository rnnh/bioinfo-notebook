---
layout: default
title: Bowtie
parent: 2. Program guides
---

# Bowtie

`bowtie` can be used to:
 - index reference FASTA nucleotide genomes/sequences
 - align FASTQ sequencing reads to those genomes/sequences

 If you want to align short reads (50bp or less), [bowtie is more suitable than bowtie2](bowtie2.md#differences-between-bowtie-and-bowtie2).

## Indexing a reference genome/sequence using `bowtie-build`

 Before aligning reads to a reference genome with `bowtie`, it must be indexed using `bowtie-bui
    ld`.
This command will create six files with the extensions `.1.ebwt`, `.2.ebwt`, `.3.ebwt`, `.4.ebwt`, `.rev.1.ebwt`, and `.rev.2.ebwt`.
These six files together are the index.
Once an index has been created, the original reference genome/sequence is no longer needed to align reads.
Here's an example `bowtie2-build` command:

```
$ bowtie-build reference_sequence.fasta index_name
```

In this command, the `reference_sequence.FASTA` is the nucleotide FASTA sequence we want to index, and `index_name` is the name of the index.
There will be six files beginning with the `index_name` in the output directory: `index_name.1.ebwt`, `index_name.2.ebwt`, `index_name.3.ebwt`, `index_name.4.ebwt`, `index_name.rev.1.ebwt`, and `index_name.rev.2.ebwt`.
There's no need to specify any of these files individually in subsequent `bowtie` commands, the `index_name` alone is enough to refer to the entire index.

## Aligning reads to an indexed genome/sequence using `bowtie`

Now that the genome has been indexed, FASTQ sequencing reads can be aligned to it.
This is done using the `bowtie` command.
Here is an example `bowtie2` command:

```
$ bowtie --no-unal --threads n --sam index_name -1 reads_1.fastq -2 reads_2.fastq output.sam
```

In this command...
 
1. **`--no-unal`** is an optional argument, meaning reads that do not align to the reference genome will not be written to `sam` output
2. **`--threads`** is the number (*n*) of processors/threads used
3. **`--sam`** specifies that the output should be written in the [SAM format](file_formats.md#sam)
4. **`index_name`** is the name of the genome index
4. **`-1`** is the file(s) containing mate 1 reads ([`reads_1.fastq`](file_formats.md#fastq))
5. **`-2`** is the file(s) containing mate 2 reads ([`reads_2.fastq`](file_formats.md#fastq))
6. **`output.sam`** is the output alignment in `sam` format

## Demonstration

In this video, `bowtie-build` is used to index `S_cere_GCF_000146045.2_R64_genomic.fna`, which is a copy of the [*Saccharomyces cerevisiae* S288C genome from RefSeq](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2).
The `bowtie` command is then used to align [*Saccharomyces cerevisiae* RNAseq reads](https://www.ncbi.nlm.nih.gov/sra/SRR11462797) to this bowtie index.

[![asciicast](https://asciinema.org/a/316272.svg)](https://asciinema.org/a/316272?autoplay=1)

## Further reading

1. The `bowtie` manual: <http://bowtie-bio.sourceforge.net/manual.shtml>
