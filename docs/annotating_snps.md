---
layout: default
title: Annotating SNPs
parent: 3. Scripts
---

# Annotating SNPs

[annotating_snps.R](../scripts/annotating_snps.R) is an `R` script that cross-references annotations of genome assemblies with VCF files containing SNPs of sequencing reads aligned against
those genome assemblies.
If a SNP falls within- or upstream of- an annotated genome feature (start codon, stop codon, CDS, etc.), the script will return that feature along with the SNP.
For this script to work, these files need to use the same sequence names: e.g. if the first sequence in the VCF is called "chrI", there should be a corresponding sequence called "chrI" in the GFF file.

## Usage

To use this script, variables need to be defined on lines 28 to 32 of the script:

- The GFF file name should be assigned to the variable `GFF_file`.
- The VCF file name should be assigned to the variable `VCF_file`.
- The VCF and GFF files should be in the directory `~/bioinfo-notebook/data/`.
- The number of lines in the VCF file header should be specified in the `VCF_header.int` variable. This is the number of lines that begin with `#` in the VCF file.
- The variable `upstream.int` is used to determine how far upstream from an annotated feature a SNP can be. This can be set to 0 if you do not want upstream SNPs to be considered. Setting it to 1000 will mean that SNPs up to 1,000 bases/1kb upstream from a feature will be annotated.
- The variable 'output_name' is used to specify the name of the output file, which should end in '.tsv' as it will be a tab-separated values text file.

## See also

- [annotating_snps.R on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/annotating_snps.R)
- [File formats used in bioinformatics](file_formats.md)
- [snp_calling.sh](snp_calling.md), a script for generating VCF files of SNPs.
- [genome_annotation_SwissProt_CDS.sh](genome_annotation_SwissProt_CDS.md), a script for generating genome annotation GFF files.
