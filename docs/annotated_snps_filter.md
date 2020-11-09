---
layout: default
title: Annotated SNPs filter
parent: 3. Scripts
---

# Annotated SNPs filter

[annotated_snps_filter.R](../scripts/annotated_snps_filter.R) is an `R` script cross-references annotated SNP files created using [annotating_snps.R](annotating_snps.md).
It takes two files created using this script, and returns unique SNPs for each file.
If a SNP in File 1 is not found at the same position on the same sequence as File 2, it is returned as a unique SNP, and vice versa.
These unique SNPs are then written to new `.tsv` files.

## Usage

To use this script, variables need to be defined on lines 21 and 22 of the script:

- Assign the name of the first annotated SNP file to be filtered to 'annotated_SNP_file_1'.
- Assign the name of the second annotated SNP file to be filtered to 'annotated_SNP_file_2'.
- These files should be in the `~/bioinfo-notebook/data/` directory.
- Optional: the name of the output files can be assigned on lines 109 and 115 respectively.

## See also

- [annotated_snps_filter.R on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/annotated_snps_filter.R)
- [annotating_snps.R](annotating_snps.md)
- [File formats used in bioinformatics](file_formats.md)
- [snp_calling.sh](snp_calling.md), a script for generating VCF files of SNPs.
- [genome_annotation_SwissProt_CDS.sh](genome_annotation_SwissProt_CDS.md), a script for generating genome annotation GFF files.
