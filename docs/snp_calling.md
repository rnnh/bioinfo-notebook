---
layout: default
title: SNP calling script
parent: 3. Scripts
---

# SNP calling script

[snp_calling.sh](../scripts/snp_calling.sh) is a `bash` shell script that downloads [FASTQ](file_formats.md) sequencing reads using [fastq-dump](fastq-dump.md), aligns them to a genome using [bowtie2](bowtie2.md), and writes variants (SNPs and indels) to a variant call format (VCF) file.

## Usage

```
snp_calling.sh [-h|--help] [-l|--log -p|--processors n] 
 
 This script downloads FASTQ sequencing reads, aligns them to a reference 
 genome, and finds genetic variants (SNPs/indels) based on this alignment, 
 which are written to a variant call format (VCF) file. This script should 
 be called from the 'bioinfo-notebook/' directory. 
 
 arguments: 
 	 -h | --help		 show this help text and exit 
 	 -l | --log		 redirect terminal output to a log file 
 	 -p | --processors	 optional: set the number (n) of processors to 
 				 use (default: 1) 
```

## See also

- [snp_calling.sh on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/snp_calling.sh)
- [File formats used in bioinformatics](file_formats.md)
- [samtools](samtools.md)
- [fastq-dump](fastq-dump.md)
- [bowtie2](bowtie2.md)
