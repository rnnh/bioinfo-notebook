#!/usr/bin/env python3
# https://github.com/rnnh/bioinfo-notebook.git
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
import argparse
import sys
import os

# Parsing command line arguments
parser = argparse.ArgumentParser(
    description = "Combines the featureCounts output tables in the target \
    directory.")

# -d PATH -o CUSTOM_FILENAME
parser.add_argument("-d", "--directory", dest = "path",
                    help = "path to target directory. \
                            Default: current directory")
parser.add_argument("-o", "--output", dest ="custom_filename",
                    help = "output filename.\
                            Default: featCounts_{species}_{date}.csv")

args = parser.parse_args()

# Changing to the target directory
if args.path is not None:
    path = args.path
else:
    path = os.getcwd()
os.chdir(path)

# Creating variables
fixed_headers = ["Geneid", "Chromosome", "Start", "End", "Strand", "Length"]
target_file_prefix = "feature_counts_"
date = strftime("%Y%m%d", gmtime())
counts_table = pd.DataFrame()
output_filename = str()
target_file_count = 0
species_name = str()
srr = str()

# Iterating through files in target directory, combining feature counts
# into one DataFrame object ("counts_table")
for filename in os.listdir():
    if filename.startswith(target_file_prefix):
        target_file_count = target_file_count + 1
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
        # Add the gene IDs and counts to the counts_table DataFrame as columns
        # if it's empty; otherwise add the counts only
        if counts_table.empty:
            counts_table = pd.concat([gene_ids, counts], axis = 1,
                                     sort = False)
        else:
            counts_table = pd.concat([counts_table, counts], axis = 1,
                                     sort = False)
        del featCounts_headers

if target_file_count == 0:
    # Exiting script if there are no target files in the target directory
    print("ERROR: There are no featureCount files in the target directory. \n")
    parser.print_help(sys.stderr)
    exit
else:
    # Exporting counts_table DataFrame as a CSV file
    if args.custom_filename is not None:
        output_filename = args.custom_filename
    else:
        output_filename = "featCounts_" + species_name + "_" + date + ".csv"
    counts_table.to_csv(output_filename, index = False)
