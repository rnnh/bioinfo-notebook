---
layout: default
title: fastq-dump
parent: 2. Program guides
nav_order: 3
---

# fastq-dump

`fastq-dump` is a tool for downloading sequencing reads from [NCBI's Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).
These sequence reads will be downloaded as FASTQ files.
How these FASTQ files are formatted depends on the `fastq-dump` options used.

## Downloading reads from the SRA using `fastq-dump`

In this example, we want to download FASTQ reads for a mate-pair library.

```
$ fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir path/to/reads/ SRR_ID
```

In this command...

1. **`--gzip`**: Compress output using gzip. Gzip archived reads can be read directly by [bowtie2](bowtie2.md).
2. **`--skip-technical`**: Dump only biological reads, skip the technical reads.
3. **`--readids`** or **`-I`**: Append read ID after spot ID as 'accession.spot.readid'. With this flag, one sequence gets appended the ID `.1` and the other `.2`. Without this option, pair-ended reads will have identical IDs.
4. **`--read-filter pass`**: Only returns reads that pass filtering (without `N`s).
5. **--dumpbase`** or **`-B`**: Formats sequence using base space (default for other than SOLiD). Included to avoid colourspace (in which pairs of bases are represented by numbers).
6. **`--split-3`** separates the read into left and right ends. If there is a left end without a matching right end, or a right end without a matching left end, they will be put in a single file.
7. **`--clip`** or **`-W`**: Some of the sequences in the SRA contain tags that need to be removed. This will remove those sequences.
8. **`--outdir`** or **`-O`**: *(Optional)* Output directory, default is current working directory.
9. **`SRR_ID`**: This is is the ID of the run from SRA to be downloaded. This ID begins with "SRR" and is followed by around seven digits (e.g. `SRA1234567`).

Other options that can be used instead of `--split-3`:

1. **`--split-files`** splits the FASTQ reads into two files: one file for mate 1s (`...1`), and another for mate 2s (`..._2`). This option will not mateless pairs into a third file.
2. **`--split-spot`** splits the FASTQ reads into two (mate 1s and mate 2s) within one file. `--split-spot` gives you an 8-line fastq format where forward precedes reverse (see <https://www.biostars.org/p/178586/#258378>).

## Demonstration

In this demo, `fastq-dump` is used to download compressed FASTQ reads.

[![asciicast](https://asciinema.org/a/306937.svg)](https://asciinema.org/a/306937?autoplay=1)

## Further reading

1. Rob Edward's notes on `fastq-dump`: <https://edwards.sdsu.edu/research/fastq-dump/>
