---
layout: default
title: SNP calling script
parent: 3. Scripts
---

# SNP calling script

[snp_calling.sh](../scripts/snp_calling.sh) is a `bash` shell script that downloads [FASTQ](file_formats.md) sequencing reads using [fastq-dump](fastq-dump.md), aligns them to a genome using [bowtie2](bowtie2.md), and writes variants (SNPs and indels) to a variant call format (VCF) file.

## Usage

```
snp_calling.sh [-h|--help] [-1|--one -2|--two -r|--reference] 
 [-d|--demo] [-o|--output -l|--log -p|--processors n] 
 
 This script aligns sequencing reads to a reference genome, and finds genetic 
 variants (SNPs/indels) based on this alignment, which are written to a variant
 call format (VCF) file.
 
 Calling this script with the argument '-d' or '--demo' will run this script 
 using Saccharomyces cerevisiae FASTQ sequencing reads and a Saccharomyces 
 cerevisiae reference genome, which will be downloaded from NCBI. 
 
 This script should be called from the 'bioinfo-notebook/' directory.The 
 programs required for this script are in the 'bioinfo-notebook' conda 
 environment (bioinfo-notebook/envs/bioinfo-notebook.yml or 
 bioinfo-notebook/envs/bioinfo-notebook.txt). 
 If the input files are not in the 'bioinfo-notebook/data/' directory, the full 
 file paths should be given.

 
 arguments: 
 	 -h | --help		 show this help text and exit 
 	 -1 | --one		 forward reads to align with reference sequence 
 				 (FASTQ: .fastq or .fastq.gz) 
 	 -2 | --two		 reverse reads to align with reference sequence 
 				 (FASTQ: .fastq or .fastq.gz) 
 	 -r | --reference	 reference sequence to align reads against 
 				 (FASTA nucleotide file: .fna) 
 	 -d | --demo		 run the script with demonstration inputs
 
 optional arguments: 
 	 -o | --output		 optional: name of final output file 
 				 (default: 'reference_seq_vs_reads_var.vcf', or 
 				 'S_cere_DRR237290_var.vcf' if demo is used). 
 	 -l | --log		 redirect terminal output to a log file in the 
 				 directory bioinfo-notebook/results/ 
 	 -p | --processors	 optional: set the number (n) of processors to 
 				 use (default: 1)
```

## See also

- [snp_calling.sh on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/snp_calling.sh)
- [File formats used in bioinformatics](file_formats.md)
- [samtools](samtools.md)
- [fastq-dump](fastq-dump.md)
- [bowtie2](bowtie2.md)
