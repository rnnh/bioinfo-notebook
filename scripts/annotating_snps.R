# https://github.com/rnnh/bioinfo-notebook.git

# Aim ==========================================================================

# The aim of this script is to cross-reference annotations of genome assemblies
# with VCF files containing SNPs of sequencing reads aligned against those
# genome assemblies. If a SNP falls within- or upstream of- an annotated
# genome feature (start codon, stop codon, CDS, etc.), the script will return
# that feature along with the SNP.

# Selecting files and parameters ===============================================

#   - The VCF and GFF files to be cross-referenced are specified in this
#     section. For this script to work, these files need to use the same
#     sequence names: e.g. if the first sequence in the VCF is called "chrI",
#     there should be a corresponding sequence called "chrI" in the GFF file.
#   - The VCF and GFF files should be in the directory '/bioinfo-notebook/data/'.
#   - The number of lines in the VCF file header should be specified in the
#     'VCF_header.int' variable. This is the number of lines that begin with '#'
#     in the VCF file.
#   - The variable 'upstream.int' is used to determine how far upstream from an
#     annotated feature a SNP can be. This can be set to 0 if you do not want
#     upstream SNPs to be considered. Setting it to 1000 will mean that SNPs
#     up to 1,000 bases/1kb upstream from a feature will be annotated.
#   - The variable 'output_name' is used to specify the name of the output file,
#     which should end in '.tsv' as it will be a tab-separated values text file.

GFF_file <- "<.gff File name here>"
VCF_file <- "<.vcf File name here>"
VCF_header.int <- as.integer("<Number of header lines in .vcf file here>")
upstream.int <- as.integer("<Number of bases upstream of a feature a SNP can be")
output_name <- "Output_file_name_here.tsv"

# Setup ========================================================================

# Loading dplyr
if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")
library(dplyr)

# Setting the working directory
setwd("~/bioinfo-notebook/data")

# Defining headers for GFF and VCF files
gff_headers <- c("sequence", "source", "feature", "start", "end", "score",
                 "strand", "phase", "attributes")
vcf_headers <- c("sequence", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                 "FORMAT", "SAMPLE")

# Reading in files and applying format headers =================================

# Reading assembly annotation file (.gff)
genome_annotation.df <- read.delim(
  GFF_file,
  header = FALSE, stringsAsFactors = FALSE
)
colnames(genome_annotation.df) <- gff_headers

# SNPs from alignment with annotated genome (.vcf)
SNPs.df <- read.delim(
  VCF_file,
  header = FALSE, stringsAsFactors = FALSE, skip = VCF_header.int
)
colnames(SNPs.df) <- vcf_headers

# Joining data frames using dplyr ==============================================

# In this section...
#   - inner_join() is used to join together the genome annotation and SNPs
#     data frames along the column "sequence": i.e. rows with the same value
#     for "sequence" are joined together
#   - select() is used to remove the columns ID, FORMAT and FILTER
#   - filter() is used to filter out SNPs which do not fall within regions of
#     interest: i.e. SNPs that are not within- or upstream of- features

# Joining data frames with genome annotation and SNPs
SNPs_with_annotations.df <- inner_join(genome_annotation.df, SNPs.df, 
                                       by = "sequence") %>%
                                       select(-ID, -FORMAT, -FILTER) %>%
                                       filter(POS >= (start - upstream) &
                                              POS <= end)
                                              
# Removing redundant data frames
rm(genome_annotation.df, SNPs.df)

# Ordering filtered data frame of SNPs with annotations ========================
attach(SNPs_with_annotations.df)
SNPs_with_annotations.df <- SNPs_with_annotations.df[order(sequence, start, end), ]
detach(SNPs_with_annotations.df)

# Exporting SNPs with annotations to tab-separated value (.tsv) file ===========
write.table(SNPs_with_annotations.df,
            file = "DMKU3_annotation_with_ATCC_SNPs.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)

# Exiting ======================================================================
quit(save = "no")
