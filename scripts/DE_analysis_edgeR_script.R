# Loading required libraries
library(tidyverse)
library(limma)
library(edgeR)

# Changing working directory
# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

# Reading in the feature count file as "counts.df"
counts.df <- read.csv("../data/featCounts_S_cere_20200331.csv")

# Printing the start of the counts.df object in R...
head(counts.df)

# Using the "Geneid" column to set the rownames
rownames(counts.df) <- counts.df$Geneid

# Removing the "Geneid" column
counts.df$Geneid <- NULL

# Printing the start of the counts.df object in R...
head(counts.df)

# Reading in the design table as "design.df"
design.df <- read.csv("../data/design_table.csv", fileEncoding="UTF-8-BOM")

# Printing the start of the design.df object in R...
print(design.df)

# Subsetting gene counts according to experimental condition
counts_standard.df  <- counts.df[,c("SRR8933535", "SRR8933536", "SRR8933537")]
counts_anaerobic.df <- counts.df[,c("SRR8933506", "SRR8933511", "SRR8933512")]
counts_high_temp.df <- counts.df[,c("SRR8933532", "SRR8933533", "SRR8933534")]
counts_low_pH.df    <- counts.df[,c("SRR8933530", "SRR8933531", "SRR8933539")]
counts_pressure.df  <- counts.df[,c("SRR8933509", "SRR8933510", "SRR8933538")]

# Printing the structure of the gene counts set and subsets
str(counts.df)
str(counts_standard.df)
str(counts_anaerobic.df)
str(counts_high_temp.df)
str(counts_low_pH.df)
str(counts_pressure.df)

# Defining function "RSD.test()"
RSD.test <- function(dataframe){
  # This function tests whether the relative standard deviation (RSD) is less
  # than or equal to one for each row in a data frame.
  # It adds the result to a new variable in the data frame called "RSD.test".
  # For a given row, if data.frame$RSD.test is TRUE, that row has an RSD less
  # than or equal to one, i.e. RSD <= 1.
  # If data.frame$RSD.test is FALSE, that row has an RSD outside of this range.
  RSD_tests = dataframe[,1]
  for (row_index in 1:nrow(dataframe)){
    row = as.numeric(dataframe[row_index,])
    RSD = sd(row) / mean(row)
    RSD_tests[row_index] = RSD <= 1 || is.na(RSD)
  }
  dataframe$RSD.test <- as.factor(RSD_tests)
  levels(dataframe$RSD.test) <- c(FALSE, TRUE)
  return(dataframe)
}

# Applying RSD.test() to gene count subsets
counts_standard.df  <- RSD.test(counts_standard.df)
counts_anaerobic.df <- RSD.test(counts_anaerobic.df)
counts_high_temp.df <- RSD.test(counts_high_temp.df)
counts_low_pH.df    <- RSD.test(counts_low_pH.df)
counts_pressure.df  <- RSD.test(counts_pressure.df)

# Printing the structure of the gene counts subsets
str(counts_standard.df)
str(counts_anaerobic.df)
str(counts_high_temp.df)
str(counts_low_pH.df)
str(counts_pressure.df)

# Creating list of genes which failed RSD test
RSD_failed_genes <- rownames(counts_standard.df[
  which(counts_standard.df$RSD.test == FALSE),])
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_anaerobic.df[
  which(counts_anaerobic.df$RSD.test == FALSE),]))
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_high_temp.df[
  which(counts_high_temp.df$RSD.test == FALSE),]))
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_low_pH.df[
  which(counts_low_pH.df$RSD.test == FALSE),]))
RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_pressure.df[
  which(counts_pressure.df$RSD.test == FALSE),]))
RSD_failed_genes <- unique(RSD_failed_genes)
length(RSD_failed_genes)

# Filtering gene counts
filtered_counts.df <- counts.df[
  which(!rownames(counts.df) %in% RSD_failed_genes),]

# Printing the structure of the filtered gene counts
str(filtered_counts.df)

# Checking that gene counts were correctly filtered
nrow(counts.df) - length(RSD_failed_genes) == nrow(filtered_counts.df)

# Removing redundant objects from R environment
rm(counts_anaerobic.df, counts_high_temp.df, counts_low_pH.df,
   counts_pressure.df, counts_standard.df, counts.df, RSD_failed_genes)

# Creating a DGEList object using the filtered gene counts
counts.DGEList <- DGEList(counts = filtered_counts.df,
                          genes = rownames(filtered_counts.df))

# Printing the design table
print(design.df)

# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df) == design.df$run)

# Add grouping information to DGEList object
counts.DGEList$samples$group <- as.factor(design.df$condition)

# Printing counts.DGEList
counts.DGEList

# Getting counts per million (CPM) for each gene
counts.cpm  <- cpm(counts.DGEList)

# Summary of the counts.DGEList object: number of genes, number of samples
dim(counts.DGEList)

# Creating an object to filter genes with a low number of reads
low_read_filter <- rowSums(counts.cpm) >= 2
summary(low_read_filter)

# Filtering genes with a low number of reads
counts.DGEList <- counts.DGEList[rowSums(counts.cpm) >= 2, ]
dim(counts.DGEList)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in low_read_filter
length(low_read_filter[low_read_filter == TRUE]) == dim(counts.DGEList)[1]

# Removing the low read filter
rm(low_read_filter)

# Normalising counts
# Printing library size per sample
counts.DGEList$samples$lib.size
# Printing the normalisation factors for the libraries
counts.DGEList$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
counts.DGEList <- calcNormFactors(counts.DGEList, method = "TMM")
counts.DGEList$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition_ <- design.df$condition
counts.DGEList <- estimateDisp(counts.DGEList,
                               design = model.matrix(~condition_))

counts.DGEList

# Exact tests for differences between experimental conditions
condition_

std_anaerobic.DGEExact <- exactTest(counts.DGEList, pair = c("standard",
                                                             "anaerobic"))
std_salt.DGEExact <- exactTest(counts.DGEList, pair = c("standard",
                                                        "osmotic_pressure"))
std_temp.DGEExact <- exactTest(counts.DGEList, pair = c("standard",
                                                        "high_temp"))
std_pH.DGEExact <- exactTest(counts.DGEList, pair = c("standard",
                                                      "low_pH"))

# Extracting most differentially expressed genes from exact tests
std_anaerobic.topTags <- topTags(std_anaerobic.DGEExact)
std_salt.topTags <- topTags(std_salt.DGEExact)
std_temp.topTags <- topTags(std_temp.DGEExact)
std_pH.topTags <- topTags(std_pH.DGEExact)

# Printing the most differentially expressed genes
std_anaerobic.topTags
std_salt.topTags
std_temp.topTags
std_pH.topTags

# Printing session information
sessionInfo()