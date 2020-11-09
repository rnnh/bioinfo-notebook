---
layout: default
title: SPAdes
parent: 2. Program guides
---

# SPAdes

SPAdes [is an assembly toolkit containing various assembly pipelines](https://github.com/ablab/spades/blob/spades_3.14.1/README.md).

## Assembling a genome from Illumina paired-end reads using `SPAdes`

`SPAdes` can be used to assemble paired-end reads as follows:

```bash
$ spades -1 reads_1.fq.gz -2 reads_2.fq.gz -t 5 -m 200 -o results/directory/
```

In this command...

1. **`-1`** is the file with forward reads.
2. **`-2`** is the file with reverse reads.
3. **`-t`** or **`--threads`** sets the number of processors/threads to use. The default is 16.
4. **`-m`** or **`--memory`** is memory the limit in Gb. SPAdes terminates if it reaches this limit. The default value is 250Gb.
5. **`-o`** or **`--outdir`** is the output directory to use. The default is the current directory.

SPAdes supports uncompressed (**`.fastq`** or **`.fq`**) or compressed (**`.fastq.gz`** or **`.fq.gz`**) sequencing read inputs.
In the output directory, the assembled genome will be available as contigs (**`contigs.fasta`**) and scaffolds (**`scaffolds.fasta`**), both of which are FASTA nucleotide files.

## See also

- [conda](conda.md): The `bioinfo-notebook` conda environment includes SPAdes
- [File formats used in bioinformatics](file_formats.md)
- [snp_calling.sh](snp_calling.md): This script uses SPAdes

## Further reading

- [SPAdes on GitHub](https://github.com/ablab/spades/)
