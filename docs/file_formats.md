# File formats

*work in progress*

A brief introduction to various file formats used in bioinformatics.

## Table of contents

- [Sequence formats](#sequence-formats)
    - [FASTA](#fasta)
    - [FASTQ](#fastq)
- [Alignment formats](#alignment-formats)
    - [BAM](#bam)
    - [CRAM](#cram)
    - [SAM](#sam)
- [Genome annotation formats](#genome-annotation-formats)
    - [GFF](#gff)
    - [GTF](#gtf)

## Sequence formats

These are file formats for storing nucleotide sequences and/or amino acid (protein) sequences.

### FASTA

FASTA is a ubiquitous text-based format for representing nucleotide sequences or amino acid sequences. A FASTA file can contain one sequence or multiple sequences. If a FASTA file contains multiple sequences, it may sometimes be referred to as a "multi-FASTA" file.

Here is an example of a FASTA file:

```
>gi|1817694395|ref|NZ_JAAGMU010000151.1| Streptomyces sp. SID7958 contig-52000002, whole genome shotgun sequence
CCGGCTGGCGCGGCTGGCGCTGGCGGTGGGGCTGCGGCTGCTGGAGCTGGGGGTGGCGCTGGAGGCGCAC
GGCCAGAACCTGCTGGTGGTGCTGTCGCCGTCCGGGGAGCCGCGGCGGCTGGTCTACCGCGATCTGGCGG
ACATCCGGGTCTCCCCCGCGCGGCTGGCCCGGCACGGTATCCGGGTTCCGGACCTGCCGGCG

>gi|1643051563|gb|SZWM01000399.1| Citrobacter sp. TBCS-14 contig3128, whole genome shotgun sequence
GCACAGTGAGATCAGCATTCCGTTGGATCTACTGGTCAATCAAAACCTGACGCTGGGTACTGAATGGAAC
CAGCAGCGCATGAAGGACATGCTGTCTAACTCGCAGACCTTTATGGGCGGTAATATTCCAGGCTACAGCA
GCACCGATCGCAGCCCATATTCGAAAGCCGAGATCTTCTCTTTGTTTGCCGAAAACAACATG
```

This example FASTA file contains two linear nucleotide sequences. Each FASTA entry begins with a `>` (greater-than) symbol, followed by a comment on the same line describing the sequence that will follow. The actual sequence begins on the line after this comment. Another `>` (greater-than) symbol denotes the beginning of another FASTA entry (comment describing the sequence + the sequence itself).

#### FASTA filename extensions

FASTA files usually end with the extension `.fasta`. This extension is arbitrary, as the content of the file determines its format, not its extension. More descriptive filename extensions can be used instead of `.fasta`, which are useful as they describe the type of sequence(s) in the file at a glance. 

Here are some examples...
- **`.fna`** can be used for **F**ASTA **n**ucleic **a**cids
- **`.faa`** can be used for **F**ASTA **a**mino **a**cids
- **`.frn`** can be used for **F**ASTA non-coding **RN**A

### FASTQ

The FASTQ format is an extension of [FASTA](#fasta) that stores both biological sequences (usually nucleotide sequences) and their corresponding quality scores. Both the sequence letter and quality score are encoded with a single character for brevity.

A FASTQ file normally uses four lines per sequence:
1. A line beginning with `@` followed by a sequence identifier and optional description (like the comment line at in a FASTA file)
2. The raw sequence letters
3. A line beginning with `+`, sometimes followed by the same comment as the first line
4. A line encoding the quality values for the sequence in line 2, with the same numbers of symbols as letters in the sequence

Here's an example FASTQ file:

```
@SRR8933535.1 1 length=75
NAGGAAACAAAGGCTTACCCGTTATCATTTCCGCAAGAATGCACCCACACGACCATATATCAATGGATGTGGAGT
+SRR8933535.1 1 length=75
#AAAAEEEEEEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAE
@SRR8933535.2 2 length=75
NGAGGAGTGGTGGTAGTGTTGCTTGGTGGCAAAGATGTAGTTGGTGGGAAAGCTGAAGTGGTACCGTTGGTTGGA
+SRR8933535.2 2 length=75
#AAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEE
```
