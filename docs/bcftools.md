---
layout: default
title: Bcftools
parent: 2. Program guides
---

# Bcftools

Bcftools are a set of [utilities for variant calling and manipulating VCFs and BCFs](https://samtools.github.io/bcftools/bcftools.html).

## Generating genotype likelihoods for alignment files using `bcftools mpileup`

`bcftools mpileup` can be used to generate VCF or BCF files containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files as follows:

```bash
$ bcftools mpileup --max-depth 10000 --threads n -f reference.fasta -o genotype_likelihoods.bcf reference_sequence_alignmnet.bam
```

In this command...

1. **`--max-depth`** or **`-d`** sets the reads per input file for each position in the alignment. In this case, it is set to 10000
2. **`--threads`** sets the number (*n*) of processors/threads to use.
3. **`--fasta-ref`** or **`-f`** is used to select the [faidx-indexed FASTA](samtools.md#indexing-a-fasta-file-using-samtools-faidx) nucleotide reference file (*reference.fasta*) used for the alignment.
4. **`--output `** or **`-o`** is used to name the ouput file (*genotype_likelihoods.bcf*).
5. The final argument given is the input BAM alignment file (*reference_sequence_alignment.bam*). Multiple input files can be given here.

## Variant calling using `bcftools call`

`bcftools call` can be used to call SNP/indel variants from a BCF file as follows:

```bash
$ bcftools call -O b --threads n -vc --ploidy 1 -p 0.05 -o variants_unfiltered.bcf genotype_likelihoods.bcf
```

In this command...

1. **`--output-type`** or **`-O`** is used to select the output format. In this case, *b* for BCF.
2. **`--threads`** sets the number (*n*) of processors/threads to use.
3. **`-vc`** specifies that we want the output to contain variants only, using the original [SAMtools](samtools.md) consensus caller.
4. **`--ploidy`** specifies the ploidy of the assembly.
5. **`--pval-threshold`** or **`-p`** is used to the set the p-value threshold for variant sites (*0.05*).
6. **`--output `** or **`-o`** is used to name the ouput file (*variants_unfiltered.bcf*).
7. The final argument is the input BCF file (*genotype_likelihoods.bcf*).

## Filtering variants using `bcftools filter`

`bcftools filter` can be used to filter variants from a BCF file as follows...

```bash
$ bcftools filter --threads n -i '%QUAL>=20' -O v -o variants_filtered.vcf variants_unfiltered.bcf
```

In this command...

1. **`--threads`** sets the number (*n*) of processors/threads to use.
2. **`--include`** or **`-i`** is used to define the expression used to filter sites. In this case, *`%QUAL>=20`* results in sites with a quality score greater than or equal to 20.
3. **`--output-type`** or **`-O`** is used to select the output format. In this case, *v* for VCF.
4. **`--output `** or **`-o`** is used to name the ouput file (*variants_filtered.vcf*).
5. The final argument is the input BCF file (*genotype_likelihoods.bcf*).

## See also

- [File formats used in bioinformatics](file_formats.md)
- [SNP calling script](snp_calling.md)

## Futher reading

- [bcftools documentation](https://samtools.github.io/bcftools/bcftools.html)
