# File formats

A brief introduction to various file formats used in bioinformatics.

## Table of contents

- [Sequence formats](#sequence-formats)
    - [FASTA](#fasta)
    - [FASTQ](#fastq)
- [Alignment formats](#alignment-formats)
    - [SAM](#sam)
    - [BAM](#bam)
    - [CRAM](#cram)
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

## Alignment formats

These are formats for storing alignments of nucleotide or amino acid sequences.

### SAM

**S**equence **A**lignment/**M**ap (**SAM**) is a text-based alignment format that supports single- and paired-end reads produced by different sequencing platforms.
It can support short and long reads (up to 128Mbp).
The format has been extended to include unmapped sequences, and it may contain other data such as base-call & alignment qualities.

The SAM format consists of a header and an alignment section.
Headings begin with the `@` symbol, which distinguishes them from the alignment section.
All lines in a SAM file are tab-delimited.
Alignment sections contain 11 mandatory fields, with other fields being optional.
Although the mandatory fields must be present, their value can be a `*` or a `0` depending on the field.
Optional fields are presented as key-value pairs in the format `TAG:TYPE:VALUE`.

The mandatory fields in a SAM file are...

1. **`QNAME`** **Q**uery template **name**. Reads/segments with the same QNAME are from the same template. A read may occupy multiple alignment lines when its alignment is chimeric or when multiple mappings are given.
2. **`FLAG`** Bitwise **flag** (pairing, strand, mate strand, etc.).
3. **`RNAME`** **R**eference sequence **name**.
4. **`POS`** 1-based leftmost mapping **pos**ition. The first base in a reference sequence has coordinate 1. `POS` is set to 0 for an unmapped read.
5. **`MAPQ`** **Map**ping **q**uality. If equal to 255, the mapping quality is not available.
6. **`CIGAR`** Extended **C**onsise **I**diosyncratic **G**apped **A**lignment **R**eport string. `M` for match/mismatch, `I` for insertion and `D` for deletion compared with the reference, `N` for skipped bases on the reference, `S` for soft clipping, `H` for hard clipping, and `P` for padding.
7. **`RNEXT`** **R**eference name of the mate/**next** read in the template. For the last read, the next read is the first read in the template.
8. **`PNEXT`** **P**osition of the primary alignment of the mate/**next** read.
9. **`TLEN`** Observed **t**emplate **len**gth.
10. **`SEQ`** Segment **seq**uence.
11. **`QUAL`** Phred-based base **qual**ity (same as the quailty string in [FASTQ](#fastq)) plus 33.

SAM files can be manipulated using [SAMtools](samtools.md).

### BAM

**B**inary **A**lignment **M**ap is a binary representation of [SAM](#sam), containing the same information in binary format for improved performance.
A position-sorted BAM file can be indexed to allow random access.

## CRAM

CRAM is a sequencing read file format that is highly space efficient by using reference-based compression of sequence data and offers both lossless and lossy modes of compression.
CRAM files are typically 30 to 60% smaller than their BAM equivalents.
CRAM has the following major objectives:

1. Significantly better lossless compression than BAM
2. Full compatibility with BAM
3. Effortless transition to CRAM from using BAM files
4. Support for controlled loss of BAM data

## References

- [CRAM format specification (version 3.0)](https://github.com/samtools/hts-specs/blob/5a5d05fa157c679f34db8920ce3acab1d9f3dfd1/CRAMv3.pdf)
- [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352)
