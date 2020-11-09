---
layout: default
title: Fasterq-dump
parent: 2. Program guides
---

# Fasterq-dump

`fasterq-dump` is a tool for downloading sequencing reads from [NCBI's Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).
These sequence reads will be downloaded as [FASTQ files](file_formats.md#fastq).
`fasterq-dump` is a newer, streamlined alternative to [fastq-dump](fastq-dump.md); both of these programs are a part of [sra-tools](https://anaconda.org/bioconda/sra-tools).

## `fasterq-dump` vs `fastq-dump`

Here are a few of the differences between `fastq-dump` and `fasterq-dump`:

1. In `fastq-dump`, the flag `--split-3` is required to separate paired reads into left and right ends. This is the default setting in `fasterq-dump`.
2. The `fastq-dump` flag `--skip-technical` is no longer required to skip technical reads in `fasterq-dump`. Instead, the flag `--include-technical` is required to include technical reads when using `fasterq-dump`.
3. There is no `--gzip` or `--bzip2` flag in `fasterq-dump` to download compressed reads with `fasterq-dump`. However, FASTQ files downloaded using `fasterq-dump` can still be subsequently compressed.

The following commands are equivalent, but will be executed faster using `fasterq-dump`:

```
$ fastq-dump SRR_ID --split-3 --skip-technical
$ fasterq-dump SRR_ID
```

## Downloading reads from the SRA using `fasterq-dump`

In this example, we want to download FASTQ reads for a mate-pair library.

```
fastq-dump --threads n --progress SRR_ID
```

In this command...

1. **`--threads`** specifies the number (*`n`*) processors/threads to be used.
2. **`--progress`** is an optional argument that displays a progress bar when the reads are being downloaded.
3. **`SRR_ID`** is the ID of the run from the SRA to be downloaded. This ID begins with "SRR" and is followed by around seven digits (e.g. `SRA1234567`).

## Demonstration

In this video, `fasterq-dump` is used to download [*Saccharomyces cerevisiae* RNAseq reads](https://www.ncbi.nlm.nih.gov/sra/SRR11462797) from the SRA.

[![asciicast](https://asciinema.org/a/316273.svg)](https://asciinema.org/a/316273?autoplay=1)

## See also

- [fastq-dump](fastq-dump.md)

## References

1. [How to use fasterq-dump from the sra-tools wiki on GitHub](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
