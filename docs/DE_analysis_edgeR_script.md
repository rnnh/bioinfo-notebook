---
layout: default
title: DE analysis edgeR script.R
parent: 3. Scripts
---

# DE analysis edgeR script.R

In this script, differential expression (DE) analysis is carried out on RNA-seq data, using the `R` programming language with the `edgeR` library.

```{R}
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
```

```
##       Geneid SRR8933532 SRR8933534 SRR8933509 SRR8933530 SRR8933511 SRR8933533
## 1    YAL068C          7          6          8          7         24          4
## 2  YAL067W-A          0          0          0          0          0          0
## 3    YAL067C          0          0          0          0          0          0
## 4    YAL065C          0          0          1          1          1          1
## 5  YAL064W-B         23         28         24         19         32         24
## 6  YAL064C-A          0          0          0          0          0          0
##   SRR8933537 SRR8933506 SRR8933531 SRR8933538 SRR8933512 SRR8933510 SRR8933535
## 1          6         11         12          3         28          5          5
## 2          0          0          0          0          0          0          0
## 3          0          0          0          0          0          0          0
## 4          3          3          2          3          0          0          2
## 5         11         17         11         19         18         25         17
## 6          0          0          0          0          0          0          0
##   SRR8933536 SRR8933539
## 1          6          5
## 2          0          0
## 3          0          0
## 4          1          1
## 5         15         15
## 6          0          0
```

```{R}
# Using the "Geneid" column to set the rownames
rownames(counts.df) <- counts.df$Geneid

# Removing the "Geneid" column
counts.df$Geneid <- NULL

# Printing the start of the counts.df object in R...
head(counts.df)
```

```
##           SRR8933532 SRR8933534 SRR8933509 SRR8933530 SRR8933511 SRR8933533
## YAL068C            7          6          8          7         24          4
## YAL067W-A          0          0          0          0          0          0
## YAL067C            0          0          0          0          0          0
## YAL065C            0          0          1          1          1          1
## YAL064W-B         23         28         24         19         32         24
## YAL064C-A          0          0          0          0          0          0
##           SRR8933537 SRR8933506 SRR8933531 SRR8933538 SRR8933512 SRR8933510
## YAL068C            6         11         12          3         28          5
## YAL067W-A          0          0          0          0          0          0
## YAL067C            0          0          0          0          0          0
## YAL065C            3          3          2          3          0          0
## YAL064W-B         11         17         11         19         18         25
## YAL064C-A          0          0          0          0          0          0
##           SRR8933535 SRR8933536 SRR8933539
## YAL068C            5          6          5
## YAL067W-A          0          0          0
## YAL067C            0          0          0
## YAL065C            2          1          1
## YAL064W-B         17         15         15
## YAL064C-A          0          0          0
```

```{R}
# Reading in the design table as "design.df"
design.df <- read.csv("../data/design_table.csv", fileEncoding="UTF-8-BOM")

# Printing the start of the design.df object in R...
print(design.df)
```

```
##           run          name         condition
## 1  SRR8933532  SCEhightemp3         high_temp
## 2  SRR8933534  SCEhightemp1         high_temp
## 3  SRR8933509       SCEkcl3  osmotic_pressure
## 4  SRR8933530     SCElowPH2            low_pH
## 5  SRR8933511     SCEanaer2         anaerobic
## 6  SRR8933533  SCEhightemp2         high_temp
## 7  SRR8933537      SCEstan1          standard
## 8  SRR8933506     SCEanaer3         anaerobic
## 9  SRR8933531     SCElowPH1            low_pH
## 10 SRR8933538       SCEkcl1  osmotic_pressure
## 11 SRR8933512     SCEanaer1         anaerobic
## 12 SRR8933510       SCEkcl2  osmotic_pressure
## 13 SRR8933535      SCEstan3          standard
## 14 SRR8933536      SCEstan2          standard
## 15 SRR8933539     SCElowPH3            low_pH
```

```{R}
# Subsetting gene counts according to experimental condition
counts_standard.df  <- counts.df[,c("SRR8933535", "SRR8933536", "SRR8933537")]
counts_anaerobic.df <- counts.df[,c("SRR8933506", "SRR8933511", "SRR8933512")]
counts_high_temp.df <- counts.df[,c("SRR8933532", "SRR8933533", "SRR8933534")]
counts_low_pH.df    <- counts.df[,c("SRR8933530", "SRR8933531", "SRR8933539")]
counts_pressure.df  <- counts.df[,c("SRR8933509", "SRR8933510", "SRR8933538")]
```

```{R}
# Printing the structure of the gene counts set and subsets
str(counts.df)
```

```
## 'data.frame':	6420 obs. of  15 variables:
##  $ SRR8933532: int  7 0 0 0 23 0 0 0 26 1124 ...
##  $ SRR8933534: int  6 0 0 0 28 0 0 0 25 1045 ...
##  $ SRR8933509: int  8 0 0 1 24 0 0 0 30 556 ...
##  $ SRR8933530: int  7 0 0 1 19 0 0 0 53 1135 ...
##  $ SRR8933511: int  24 0 0 1 32 0 0 0 66 252 ...
##  $ SRR8933533: int  4 0 0 1 24 0 0 0 30 1081 ...
##  $ SRR8933537: int  6 0 0 3 11 0 0 0 33 1288 ...
##  $ SRR8933506: int  11 0 0 3 17 0 0 0 30 235 ...
##  $ SRR8933531: int  12 0 0 2 11 0 0 0 31 388 ...
##  $ SRR8933538: int  3 0 0 3 19 0 0 0 49 569 ...
##  $ SRR8933512: int  28 0 0 0 18 0 0 0 61 209 ...
##  $ SRR8933510: int  5 0 0 0 25 0 0 0 49 567 ...
##  $ SRR8933535: int  5 0 0 2 17 0 0 0 34 350 ...
##  $ SRR8933536: int  6 0 0 1 15 0 0 0 32 474 ...
##  $ SRR8933539: int  5 0 0 1 15 0 0 0 64 1236 ...
```

```{R}
str(counts_standard.df)
```

```
## 'data.frame':	6420 obs. of  3 variables:
##  $ SRR8933535: int  5 0 0 2 17 0 0 0 34 350 ...
##  $ SRR8933536: int  6 0 0 1 15 0 0 0 32 474 ...
##  $ SRR8933537: int  6 0 0 3 11 0 0 0 33 1288 ...
```

```{R}
str(counts_anaerobic.df)
```

```
## 'data.frame':	6420 obs. of  3 variables:
##  $ SRR8933506: int  11 0 0 3 17 0 0 0 30 235 ...
##  $ SRR8933511: int  24 0 0 1 32 0 0 0 66 252 ...
##  $ SRR8933512: int  28 0 0 0 18 0 0 0 61 209 ...
```

```{R}
str(counts_high_temp.df)
```

```
## 'data.frame':	6420 obs. of  3 variables:
##  $ SRR8933532: int  7 0 0 0 23 0 0 0 26 1124 ...
##  $ SRR8933533: int  4 0 0 1 24 0 0 0 30 1081 ...
##  $ SRR8933534: int  6 0 0 0 28 0 0 0 25 1045 ...
```

```{R}
str(counts_low_pH.df)
```

```
## 'data.frame':	6420 obs. of  3 variables:
##  $ SRR8933530: int  7 0 0 1 19 0 0 0 53 1135 ...
##  $ SRR8933531: int  12 0 0 2 11 0 0 0 31 388 ...
##  $ SRR8933539: int  5 0 0 1 15 0 0 0 64 1236 ...
```

```{R}
str(counts_pressure.df)
```

```
## 'data.frame':	6420 obs. of  3 variables:
##  $ SRR8933509: int  8 0 0 1 24 0 0 0 30 556 ...
##  $ SRR8933510: int  5 0 0 0 25 0 0 0 49 567 ...
##  $ SRR8933538: int  3 0 0 3 19 0 0 0 49 569 ...
```

```{R}
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
```

```{R}
# Printing the structure of the gene counts subsets
str(counts_standard.df)
```

```
## 'data.frame':	6420 obs. of  4 variables:
##  $ SRR8933535: int  5 0 0 2 17 0 0 0 34 350 ...
##  $ SRR8933536: int  6 0 0 1 15 0 0 0 32 474 ...
##  $ SRR8933537: int  6 0 0 3 11 0 0 0 33 1288 ...
##  $ RSD.test  : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 2 2 2 2 2 2 2 ...
```

```{R}
str(counts_anaerobic.df)
```

```
## 'data.frame':	6420 obs. of  4 variables:
##  $ SRR8933506: int  11 0 0 3 17 0 0 0 30 235 ...
##  $ SRR8933511: int  24 0 0 1 32 0 0 0 66 252 ...
##  $ SRR8933512: int  28 0 0 0 18 0 0 0 61 209 ...
##  $ RSD.test  : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 1 2 2 2 2 2 2 ...
```

```{R}
str(counts_high_temp.df)
```

```
## 'data.frame':	6420 obs. of  4 variables:
##  $ SRR8933532: int  7 0 0 0 23 0 0 0 26 1124 ...
##  $ SRR8933533: int  4 0 0 1 24 0 0 0 30 1081 ...
##  $ SRR8933534: int  6 0 0 0 28 0 0 0 25 1045 ...
##  $ RSD.test  : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 1 2 2 2 2 2 2 ...
```

```{R}
str(counts_low_pH.df)
```

```
## 'data.frame':	6420 obs. of  4 variables:
##  $ SRR8933530: int  7 0 0 1 19 0 0 0 53 1135 ...
##  $ SRR8933531: int  12 0 0 2 11 0 0 0 31 388 ...
##  $ SRR8933539: int  5 0 0 1 15 0 0 0 64 1236 ...
##  $ RSD.test  : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 2 2 2 2 2 2 2 ...
```

```{R}
str(counts_pressure.df)
```

```
## 'data.frame':	6420 obs. of  4 variables:
##  $ SRR8933509: int  8 0 0 1 24 0 0 0 30 556 ...
##  $ SRR8933510: int  5 0 0 0 25 0 0 0 49 567 ...
##  $ SRR8933538: int  3 0 0 3 19 0 0 0 49 569 ...
##  $ RSD.test  : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 1 2 2 2 2 2 2 ...
```

```{R}
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
```

```
## [1] 373
```

```{R}
# Filtering gene counts
filtered_counts.df <- counts.df[
  which(!rownames(counts.df) %in% RSD_failed_genes),]

# Printing the structure of the filtered gene counts
str(filtered_counts.df)
```

```
## 'data.frame':	6047 obs. of  15 variables:
##  $ SRR8933532: int  7 0 0 23 0 0 0 26 1124 1877 ...
##  $ SRR8933534: int  6 0 0 28 0 0 0 25 1045 2280 ...
##  $ SRR8933509: int  8 0 0 24 0 0 0 30 556 618 ...
##  $ SRR8933530: int  7 0 0 19 0 0 0 53 1135 2327 ...
##  $ SRR8933511: int  24 0 0 32 0 0 0 66 252 207 ...
##  $ SRR8933533: int  4 0 0 24 0 0 0 30 1081 2217 ...
##  $ SRR8933537: int  6 0 0 11 0 0 0 33 1288 1583 ...
##  $ SRR8933506: int  11 0 0 17 0 0 0 30 235 175 ...
##  $ SRR8933531: int  12 0 0 11 0 0 0 31 388 1935 ...
##  $ SRR8933538: int  3 0 0 19 0 0 0 49 569 748 ...
##  $ SRR8933512: int  28 0 0 18 0 0 0 61 209 203 ...
##  $ SRR8933510: int  5 0 0 25 0 0 0 49 567 657 ...
##  $ SRR8933535: int  5 0 0 17 0 0 0 34 350 1762 ...
##  $ SRR8933536: int  6 0 0 15 0 0 0 32 474 1647 ...
##  $ SRR8933539: int  5 0 0 15 0 0 0 64 1236 2400 ...
```

```{R}
# Checking that gene counts were correctly filtered
nrow(counts.df) - length(RSD_failed_genes) == nrow(filtered_counts.df)
```

```
## [1] TRUE
```

```{R}
# Removing redundant objects from R environment
rm(counts_anaerobic.df, counts_high_temp.df, counts_low_pH.df,
   counts_pressure.df, counts_standard.df, counts.df, RSD_failed_genes)
```

```{R}
# Creating a DGEList object using the filtered gene counts
counts.DGEList <- DGEList(counts = filtered_counts.df,
                          genes = rownames(filtered_counts.df))

# Printing the design table
print(design.df)
```

```
##           run          name         condition
## 1  SRR8933532  SCEhightemp3         high_temp
## 2  SRR8933534  SCEhightemp1         high_temp
## 3  SRR8933509       SCEkcl3  osmotic_pressure
## 4  SRR8933530     SCElowPH2            low_pH
## 5  SRR8933511     SCEanaer2         anaerobic
## 6  SRR8933533  SCEhightemp2         high_temp
## 7  SRR8933537      SCEstan1          standard
## 8  SRR8933506     SCEanaer3         anaerobic
## 9  SRR8933531     SCElowPH1            low_pH
## 10 SRR8933538       SCEkcl1  osmotic_pressure
## 11 SRR8933512     SCEanaer1         anaerobic
## 12 SRR8933510       SCEkcl2  osmotic_pressure
## 13 SRR8933535      SCEstan3          standard
## 14 SRR8933536      SCEstan2          standard
## 15 SRR8933539     SCElowPH3            low_pH
```

```{R}
# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df) == design.df$run)
```

```
##    Mode    TRUE 
## logical      15 
```

```{R}
# Add grouping information to DGEList object
counts.DGEList$samples$group <- as.factor(design.df$condition)

# Printing counts.DGEList
counts.DGEList
```

```
## An object of class "DGEList"
## $counts
##           SRR8933532 SRR8933534 SRR8933509 SRR8933530 SRR8933511 SRR8933533
## YAL068C            7          6          8          7         24          4
## YAL067W-A          0          0          0          0          0          0
## YAL067C            0          0          0          0          0          0
## YAL064W-B         23         28         24         19         32         24
## YAL064C-A          0          0          0          0          0          0
##           SRR8933537 SRR8933506 SRR8933531 SRR8933538 SRR8933512 SRR8933510
## YAL068C            6         11         12          3         28          5
## YAL067W-A          0          0          0          0          0          0
## YAL067C            0          0          0          0          0          0
## YAL064W-B         11         17         11         19         18         25
## YAL064C-A          0          0          0          0          0          0
##           SRR8933535 SRR8933536 SRR8933539
## YAL068C            5          6          5
## YAL067W-A          0          0          0
## YAL067C            0          0          0
## YAL064W-B         17         15         15
## YAL064C-A          0          0          0
## 6042 more rows ...
## 
## $samples
##                       group lib.size norm.factors
## SRR8933532        high_temp  7251197            1
## SRR8933534        high_temp  7591623            1
## SRR8933509 osmotic_pressure  7286624            1
## SRR8933530           low_pH  6690543            1
## SRR8933511        anaerobic  7552294            1
## 10 more rows ...
## 
## $genes
##               genes
## YAL068C     YAL068C
## YAL067W-A YAL067W-A
## YAL067C     YAL067C
## YAL064W-B YAL064W-B
## YAL064C-A YAL064C-A
## 6042 more rows ...
```

```{R}
# Getting counts per million (CPM) for each gene
counts.cpm  <- cpm(counts.DGEList)

# Summary of the counts.DGEList object: number of genes, number of samples
dim(counts.DGEList)
```

```
## [1] 6047   15
```

```{R}
# Creating an object to filter genes with a low number of reads
low_read_filter <- rowSums(counts.cpm) >= 2
summary(low_read_filter)
```

```
##   Mode   FALSE    TRUE 
## logical     168    5879
```

```{R}
# Filtering genes with a low number of reads
counts.DGEList <- counts.DGEList[rowSums(counts.cpm) >= 2, ]
dim(counts.DGEList)
```

```
## [1] 5879   15
```

```{R}
# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in low_read_filter
length(low_read_filter[low_read_filter == TRUE]) == dim(counts.DGEList)[1]
```

```
## [1] TRUE
```

```{R}
# Removing the low read filter
rm(low_read_filter)

# Normalising counts
# Printing library size per sample
counts.DGEList$samples$lib.size
```

```
##  [1] 7251197 7591623 7286624 6690543 7552294 7673803 7189631 7135118 7882580
## [10] 7870176 6865003 7389167 7627461 7807850 6556800
```

```{R}
# Printing the normalisation factors for the libraries
counts.DGEList$samples$norm.factors
```

```
##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```

```{R}
# Calculating normalisation factors and applying them to counts.DGEList
counts.DGEList <- calcNormFactors(counts.DGEList, method = "TMM")
counts.DGEList$samples$norm.factors
```

```
##  [1] 1.0209781 1.0926927 1.0142852 0.9584731 1.0565622 1.0837405 0.9260372
##  [8] 0.9738491 1.0012257 0.9814751 1.0499893 1.0048089 0.9216491 0.9409222
## [15] 0.9931300
```

```{R}
# Estimating common dispersion and tagwise dispersion
condition_ <- design.df$condition
counts.DGEList <- estimateDisp(counts.DGEList,
                               design = model.matrix(~condition_))

counts.DGEList
```

```
## An object of class "DGEList"
## $counts
##           SRR8933532 SRR8933534 SRR8933509 SRR8933530 SRR8933511 SRR8933533
## YAL068C            7          6          8          7         24          4
## YAL064W-B         23         28         24         19         32         24
## YAL063C           26         25         30         53         66         30
## YAL062W         1124       1045        556       1135        252       1081
## YAL061W         1877       2280        618       2327        207       2217
##           SRR8933537 SRR8933506 SRR8933531 SRR8933538 SRR8933512 SRR8933510
## YAL068C            6         11         12          3         28          5
## YAL064W-B         11         17         11         19         18         25
## YAL063C           33         30         31         49         61         49
## YAL062W         1288        235        388        569        209        567
## YAL061W         1583        175       1935        748        203        657
##           SRR8933535 SRR8933536 SRR8933539
## YAL068C            5          6          5
## YAL064W-B         17         15         15
## YAL063C           34         32         64
## YAL062W          350        474       1236
## YAL061W         1762       1647       2400
## 5874 more rows ...
## 
## $samples
##                       group lib.size norm.factors
## SRR8933532        high_temp  7251197    1.0209781
## SRR8933534        high_temp  7591623    1.0926927
## SRR8933509 osmotic_pressure  7286624    1.0142852
## SRR8933530           low_pH  6690543    0.9584731
## SRR8933511        anaerobic  7552294    1.0565622
## 10 more rows ...
## 
## $genes
##               genes
## YAL068C     YAL068C
## YAL064W-B YAL064W-B
## YAL063C     YAL063C
## YAL062W     YAL062W
## YAL061W     YAL061W
## 5874 more rows ...
## 
## $design
##   (Intercept) condition_high_temp condition_low_pH condition_osmotic_pressure
## 1           1                   1                0                          0
## 2           1                   1                0                          0
## 3           1                   0                0                          1
## 4           1                   0                1                          0
## 5           1                   0                0                          0
##   condition_standard
## 1                  0
## 2                  0
## 3                  0
## 4                  0
## 5                  0
## 10 more rows ...
## 
## $common.dispersion
## [1] 0.01852389
## 
## $trended.dispersion
## [1] 0.03254657 0.03074682 0.02881764 0.01407769 0.01843367
## 5874 more elements ...
## 
## $tagwise.dispersion
## [1] 0.04008573 0.01757259 0.05637753 0.15373788 0.01426743
## 5874 more elements ...
## 
## $AveLogCPM
## [1] 0.5953031 1.5665694 2.5486642 6.5981458 7.5622487
## 5874 more elements ...
## 
## $trend.method
## [1] "locfit"
## 
## $prior.df
## [1] 3.36033
## 
## $prior.n
## [1] 0.336033
## 
## $span
## [1] 0.3191663
```

```{R}
# Exact tests for differences between experimental conditions
condition_
```

```
##  [1] "high_temp"        "high_temp"        "osmotic_pressure" "low_pH"          
##  [5] "anaerobic"        "high_temp"        "standard"         "anaerobic"       
##  [9] "low_pH"           "osmotic_pressure" "anaerobic"        "osmotic_pressure"
## [13] "standard"         "standard"         "low_pH"
```

```{R}
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
```

```
## Comparison of groups:  anaerobic-standard 
##           genes     logFC    logCPM        PValue           FDR
## YNL117W YNL117W -5.995970  7.144851  0.000000e+00  0.000000e+00
## YLR413W YLR413W  4.652310  7.193190 1.036770e-292 3.047585e-289
## YLR174W YLR174W -5.944408  9.322594 3.394943e-246 6.652957e-243
## YDR046C YDR046C  6.225947  5.573382 7.647066e-240 1.123928e-236
## YOR348C YOR348C -6.751919  8.487139 3.206053e-239 3.769677e-236
## YOR011W YOR011W  3.414511  7.037274 1.078924e-209 1.057166e-206
## YPR151C YPR151C -4.783752  7.281862 2.282292e-204 1.916799e-201
## YGR067C YGR067C -4.250716  8.325035 9.445661e-192 6.941380e-189
## YNL237W YNL237W -4.161106  6.040136 9.496744e-191 6.203485e-188
## YKL217W YKL217W -4.698248 11.590578 4.050072e-175 2.381037e-172
```

```{R}
std_salt.topTags
```

```
## Comparison of groups:  osmotic_pressure-standard 
##           genes     logFC    logCPM        PValue           FDR
## YJL107C YJL107C  4.385630  6.567311 4.803825e-257 2.824169e-253
## YJL108C YJL108C  4.355598  6.144748 5.366865e-226 1.577590e-222
## IRT1       IRT1 -5.763462  6.961651 1.071199e-175 2.099193e-172
## YPR192W YPR192W -7.002595  5.511479 2.459025e-155 3.614153e-152
## YKL187C YKL187C -5.407868  8.916971 1.320162e-140 1.552246e-137
## YDL022W YDL022W  2.965100 10.007200 1.701233e-111 1.666925e-108
## YNL040W YNL040W -2.373204  6.153989  9.711834e-95  8.156553e-92
## YER062C YER062C  2.368862  7.163067  1.960021e-93  1.440370e-90
## YNL134C YNL134C -2.193439  8.279204  9.569668e-84  6.251120e-81
## YBL075C YBL075C -3.381949 10.422772  3.223047e-83  1.894829e-80
```

```{R}
std_temp.topTags
```

```
## Comparison of groups:  high_temp-standard 
##           genes      logFC    logCPM       PValue          FDR
## YBR001C YBR001C  1.2481811  7.195651 5.519364e-38 3.244834e-34
## YLL055W YLL055W -1.2984185  6.698127 4.792736e-26 1.408825e-22
## YDR248C YDR248C -1.0221538  6.304507 4.324947e-23 8.475454e-20
## YJR094C YJR094C  2.1237516  5.681368 6.723931e-21 9.882498e-18
## YLR213C YLR213C  1.1834041  4.703274 2.089727e-20 2.457101e-17
## YCL025C YCL025C  1.2949794  6.196320 9.829011e-19 9.630793e-16
## YNL134C YNL134C -0.9475343  8.279204 2.110110e-18 1.772191e-15
## YOR100C YOR100C  1.1463120  6.332454 2.457205e-18 1.805738e-15
## YHR092C YHR092C -1.9333146  6.106747 1.507036e-16 9.844295e-14
## YML008C YML008C -1.0407903 10.341063 2.288688e-16 1.345520e-13
```

```{R}
std_pH.topTags
```

```
## Comparison of groups:  low_pH-standard 
##               genes      logFC   logCPM       PValue          FDR
## YBR004C     YBR004C -0.8675698 5.499657 1.894999e-08 8.999696e-05
## RDN37-2     RDN37-2 -1.1460306 4.495122 3.061642e-08 8.999696e-05
## YMR321C     YMR321C -2.6358898 1.502487 1.867664e-07 3.660000e-04
## YHL010C     YHL010C  0.9757651 5.660352 4.212153e-07 3.668176e-04
## YGR174W-A YGR174W-A  0.8724173 4.983917 4.386824e-07 3.668176e-04
## YHL036W     YHL036W -0.7751659 6.004988 5.123259e-07 3.668176e-04
## YMR221C     YMR221C -0.8354038 6.310115 5.744536e-07 3.668176e-04
## YDL169C     YDL169C  1.1984822 6.792200 5.956031e-07 3.668176e-04
## YFL066C     YFL066C  0.7317283 6.283604 6.251769e-07 3.668176e-04
## YIL045W     YIL045W  0.7642446 7.399041 6.983764e-07 3.668176e-04
```

```{R}
# Printing session information
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17763)
## 
## Matrix products: default
## 
## locale:
##  [1] LC_COLLATE=English_Ireland.1252  LC_CTYPE=English_Ireland.1252   
##  [3] LC_MONETARY=English_Ireland.1252 LC_NUMERIC=C                    
##  [5] LC_TIME=English_Ireland.1252    
## 
## attached base packages:
##  [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] edgeR_3.30.3    limma_3.44.3    forcats_0.5.0   stringr_1.4.0  
##  [5] dplyr_1.0.0     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0    
##  [9] tibble_3.0.2    ggplot2_3.3.2   tidyverse_1.3.0
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.5       cellranger_1.1.0 pillar_1.4.6     compiler_4.0.2  
##  [5] dbplyr_1.4.4     tools_4.0.2      lattice_0.20-41  jsonlite_1.7.0  
##  [9] lubridate_1.7.9  lifecycle_0.2.0  gtable_0.3.0     pkgconfig_2.0.3 
## [13] rlang_0.4.7      reprex_0.3.0     cli_2.0.2        DBI_1.1.0       
## [17] rstudioapi_0.11  haven_2.3.1      withr_2.2.0      xml2_1.3.2      
## [21] httr_1.4.1       fs_1.4.2         generics_0.0.2   vctrs_0.3.1     
## [25] hms_0.5.3        locfit_1.5-9.4   grid_4.0.2       tidyselect_1.1.0
## [29] glue_1.4.1       R6_2.4.1         fansi_0.4.1      readxl_1.3.1    
## [33] modelr_0.1.8     blob_1.2.1       magrittr_1.5     splines_4.0.2   
## [37] backports_1.1.7  scales_1.1.1     ellipsis_0.3.1   rvest_0.3.5     
## [41] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.6    munsell_0.5.0   
## [45] broom_0.7.0      crayon_1.3.4    
```

## See also

- [DE_analysis_edgeR_script.R on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/DE_analysis_edgeR_script.R)
