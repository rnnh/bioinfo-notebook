# https://github.com/rnnh/bioinfo-notebook.git

# Aim ==========================================================================

# This script cross-references annotated SNP files created using
# annotating_snps.R. It takes two files created using this script, and returns
# unique SNPs for each file. If a SNP in File 1 is not found at the same
# position on the same sequence as File 2, it is returned as a unique SNP, and
# vice versa. These unique SNPs are then written to new .tsv files.

# Selecting files ==============================================================

#   - Assign the name of the first annotated SNP file to be filtered to
#     'annotated_SNP_file_1'
#   - Assign the name of the second annotated SNP file to be filtered to
#     'annotated_SNP_file_2'
#   - These files should be in the `~/bioinfo-notebook/data/` directory.
#   - Optional: the name of the output files can be assigned on lines 109 and
#     115 respectively.

annotated_SNP_file_1 <- "<.tsv File name here>"
annotated_SNP_file_2 <- "<.tsv File name here>"

# Setup ========================================================================

annotated_SNP_file_1 <- read.table(
  annotated_SNP_file_1,
  stringsAsFactors = FALSE, header = TRUE)

annotated_SNP_file_2 <- read.table(
  annotated_SNP_file_2,
  stringsAsFactors = FALSE, header = TRUE)

# Setting the working directory
setwd("/data/users/ronan/temp/scripts/")
  
# Finding rows in common between annotated SNP data frames =====================

# This needs to be carried out multiple times because the number of rows in
# each annotate SNP file differ. Two files may have a SNP in common, but it may
# not occur at the same row number

# Loops in this section are structured as follows:
# For every row index in a given data frame...
#   Get the row using the row index
#   If the SNP position for a given row index is in the other data frame...
#     Get the indices of the matching rows
#     For each index in the indices of matching row...
#       If the sequence names are the same for the matching rows...
#         Add that index to the matching row values
#         Keep only the unique indices

# Creating empty integer values for the matching SNPs
file_1_SNPs_common_with_file_2 <- integer()
file_2_SNPs_common_with_file_1 <- integer()

# Rows in common between file 1 and file 2
for (index in 1:nrow(annotated_SNP_file_2)){
  row = annotated_SNP_file_2[index, ]
  if (row$POS %in% annotated_SNP_file_1$POS){
    matching_row_indices = which(annotated_SNP_file_1$POS == row$POS)
    for (mr_index in matching_row_indices){
      if (annotated_SNP_file_1$sequence[mr_index] == row$sequence){
        file_1_SNPs_common_with_file_2 <- c(file_1_SNPs_common_with_file_2,
                                            mr_index)
        file_1_SNPs_common_with_file_2 <- unique(file_1_SNPs_common_with_file_2)
      }
    }
  }
}

# Rows in common between file 2 and file 1
for (index in 1:nrow(annotated_SNP_file_1)){
  row = annotated_SNP_file_1[index, ]
  if (row$POS %in% annotated_SNP_file_2$POS){
    matching_row_indices = which(annotated_SNP_file_2$POS == row$POS)
    for (mr_index in matching_row_indices){
      if (annotated_SNP_file_2$sequence[mr_index] == row$sequence){
        file_2_SNPs_common_with_file_1 <- c(file_2_SNPs_common_with_file_1,
                                            mr_index)
        file_2_SNPs_common_with_file_1 <- unique(file_2_SNPs_common_with_file_1)
      }
    }
  }
}

# Filtering SNPs in common between annotated SNP data frames ===================

# The matching row values produced by the loops in the previous section are
# used to subset each data frame: this is done by selecting non-matching rows

annotated_SNP_file_1_unique.df <- annotated_SNP_file_1[-file_1_SNPs_common_with_file_2, ]
annotated_SNP_file_2_unique.df <- annotated_SNP_file_2[-file_2_SNPs_common_with_file_1, ]

# Checking that the correct number of rows were filtered =======================

# If the correct number of rows were filtered, the following statements should
# all return TRUE

nrow(annotated_SNP_file_2) == nrow(annotated_SNP_file_2_unique.df) +
  length(file_2_SNPs_common_with_file_1)

nrow(annotated_SNP_file_1) == nrow(annotated_SNP_file_1_unique.df) + 
  length(file_1_SNPs_common_with_file_2)

# Writing data frames to tab-separated values (.tsv) files =====================

write.table(annotated_SNP_file_1_unique.df,
            file = c(annotated_SNP_file_1, "_filtered.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)
            
write.table(annotated_SNP_file_2_unique.df,
            file = c(annotated_SNP_file_2, "_filtered.tsv",
            fileEncoding = "UTF-8",
            sep = "\t",
            row.names = FALSE)

# Exiting ======================================================================
quit(save = "no")
