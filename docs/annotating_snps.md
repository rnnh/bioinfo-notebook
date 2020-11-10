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

## Annotated SNP format

The `.tsv` files created by this script have a combination of columns from the [GFF and VCF formats](file_formats.md) as follows...

1. `sequence` The name of the sequence where the feature is located.
2. `source` Keyword identifying the source of the feature, like a program (e.g. Augustus) or an organization (e.g. [SGD](https://www.yeastgenome.org/)).
3. `feature` The feature type name, like `gene` or `exon`. In a well-structured GFF file, all the children features always follow their parents in a single block (so all exons of a transcript are put after their parent `transcript` feature line and before any other parent transcript line).
4. `start` Genomic start of the feature, with a 1-base offset.
5. `end` Genomic end of the feature, with a 1-base offset.
6. `score` Numeric value that generally indicates the confidence of the source in the annotated feature. A value of `.` (a dot) is used to define a null value.
7. `strand` Single character that indicates the strand of the feature; it can assume the values of `+` (positive, or `5'->3'`), `-`, (negative, or `3'->5'`), `.` (undetermined).
8. `phase` Phase of coding sequence (CDS) features, indicating where the feature starts in relation to the reading frame. It can be either one of `0`, `1`, `2` (for CDS features) or `.` (for everything else).
9. `attributes` All the other information pertaining to this feature. The format, structure and content of this field is the one which varies the most between GFF formats.
10. `POS` The 1-based position of the variation on the given sequence.
11. `REF` The reference base (or bases in the case of an indel) at the given position on the given reference sequence.
12. `ALT` The list of alternative alleles at this position.
13. `QUAL` A quality score associated with the inference of the given alleles.
14. `FILTER` A flag indicating which of a given set of filters the variation has passed.
15. `INFO` An extensible list of key-value pairs (fields) describing the variation. Multiple fields are separated by semicolons with optional values in the format: <key>=<data>[,data].
16. `SAMPLE` For each (optional) sample described in the file, values are given for the fields listed in FORMAT. If multiple samples have been aligned to the reference sequence, each sample will have its own column.


## See also

- [annotating_snps.R on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/annotating_snps.R)
- [annotated_snps_filter.R](annotated_snps_filter.md)
- [File formats used in bioinformatics](file_formats.md)
- [snp_calling.sh](snp_calling.md), a script for generating VCF files of SNPs.
- [genome_annotation_SwissProt_CDS.sh](genome_annotation_SwissProt_CDS.md), a script for generating genome annotation GFF files.
