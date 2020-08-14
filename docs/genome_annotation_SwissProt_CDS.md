---
layout: default
title: Genome annotation script
parent: 3. Scripts
---

# Genome annotation SwissProt CDS

[genome annotation SwissProt CDS.sh](../scripts/genome_annotation_SwissProt_CDS.sh) is a bash script that annotates the coding sequences (CDS) in a given genome annotation.
It uses [BLAST](blast.md) and [MGKit](https://github.com/frubino/mgkit), which are included in the `bioinfo-notebook` [conda environment](conda.md).

## Usage

```
genome_annotation_SwissProt_CDS.sh [-h|--help] [-l|--log -p|--processors n 
 -e|--email] 
 
 This script...
 
 	 1. Downloads a Saccharomyces cerevisiae S288C genome assembly, and 
 	 the UniProtKB/Swiss-Prot amino acid sequences. 
 	 2. Creates a BLAST database from the downloaded Swiss-Prot sequences,
 	 and searches the S. cerevisiae genome against it using BLASTx with an
 	 E-value threshold of 1e-100. 
 	 3. Filters the BLASTx results, removing results with less than 90%
 	 identity.
 	 4. Creates a genome annotation GFF file from these BLASTx results.
 	 5. Adds information to the genome annotation from UniProt (protein
 	 names, KeGG ortholog information, EC numbers, etc.) 
 
 The end result ('S_cere.gff') is an annotation of the coding sequences (CDS) 
 in the S. cerevisiae genome that are described in UniProtKB/Swiss-Prot. 
 
 This script should be called from the 'bioinfo-notebook/' directory. 
 
 arguments: 
 	 -h | --help		 show this help text and exit 
 	 -l | --log		 redirect terminal output to a log file 
 	 -p | --processors	 optional: set the number (n) of processors to 
 				 use (default: 1) 
 	 -e | --email		 optional: contact email for UniProt queries
```

## See also

- [genome_annotation_SwissProt_CDS.sh on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/genome_annotation_SwissProt_CDS.sh)
- [BLAST](blast.md)
- [MGKit](https://github.com/frubino/mgkit)
- [Conda](conda.md)
