#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 12:08:41 2020

@author: ronan

This script creates a single CSV feature count table from the featureCounts
output tables in the target directory.

This combined feature count table can be used for differential expression
analysis (e.g. using DESeq2 or edgeR in R).
"""

# Loading required libraries
from time import gmtime, strftime
import pandas as pd
import os

# Changing to the target directory
path = "/home/ronan/bioinfo-notebook/data"
os.chdir(path)

# Creating variables
fixed_headers = ["Geneid", "Chromosome", "Start", "End", "Strand", "Length"]
target_file_prefix = "feature_counts_"
date = strftime("%Y%m%d", gmtime())
counts_table = pd.DataFrame()
species_name = str()
srr = str()

# Iterating through files in target directory, combining feature counts
# into one DataFrame object ("counts_table")
for filename in os.listdir():
    if filename.startswith(target_file_prefix):
        old_species_name = species_name
        filename_list = filename.split("_")
        srr = filename_list[2]
        species_name = filename_list[3] + "_" + filename_list[4]
        featCounts_df = pd.read_csv(filename, sep = "\t",
                                    lineterminator = '\n', skiprows = 1,
                                    header = 0)
        featCounts_headers = fixed_headers.copy()
        featCounts_headers += [srr]
        featCounts_df.columns = featCounts_headers
        gene_ids = featCounts_df["Geneid"]
        counts = featCounts_df[srr]
        if species_name == old_species_name:
            counts_table = pd.concat([counts_table, counts], axis = 1,
                                     sort = False)
        else:
            counts_table = pd.concat([gene_ids, counts], axis = 1,
                                     sort = False)
        del featCounts_headers
    
# Exporting counts_table DataFrame as a CSV file
output_filename = "featCounts_" + species_name + "_" + date + ".csv"
counts_table.to_csv(output_filename, index = False)
