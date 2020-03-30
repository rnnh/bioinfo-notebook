# combining featCount tables.py

This is a Python script that creates a single CSV feature count table from the featureCounts output tables in the target directory.
This combined feature count table can be used for differential expression analysis (e.g. using DESeq2 or edgeR in R).

## Demonstration

This is a video demonstartion of [combining_featCount_tables.py](scripts/combining_featCount_tables.py).

[![asciicast](https://asciinema.org/a/311771.svg)](https://asciinema.org/a/311771?autoplay=1)

In this video, `combining_featCount_tables.py` is used to combine the following [featureCounts](featureCounts.md) tables:

```
feature_counts_SRR8933506_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933509_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933510_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933511_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933512_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933530_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933531_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933532_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933533_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933534_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933535_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933536_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933537_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933538_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
feature_counts_SRR8933539_S_cere_GCF_000146045.2_R64_genomic.fna.tsv
```

These featureCounts results were generated using the following [fastq-dump_to_featureCounts.sh](fastq-dump_to_featureCounts.md) command:

```bash
$ bash ../scripts/fastq-dump_to_featureCounts.sh -a S_cere_GCF_000146045.2_R64_genomic.gtf -f S_cere_GCF_000146045.2_R64_genomic.fna --verbose -p 3 SRR8933506 SRR8933509 SRR8933510 SRR8933511 SRR8933512 SRR8933530 SRR8933531 SRR8933532 SRR8933533 SRR8933534 SRR8933535 SRR8933536 SRR8933537 SRR8933538 SRR8933539
```

In this command, the full genome sequence (`S_cere_GCF_000146045.2_R64_genomic.fna`) and genome annotation (`S_cere_GCF_000146045.2_R64_genomic.gtf`) for [*Saccharomyces cerevisiae* S288C](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2) are used.

These featureCounts results were then combined using the following command:

```bash
$ python ../scripts/combining_featCount_tables.py
```

Running this script combines all the featureCounts results in a directory into a single CSV file.
If a custom name for this file is not given, it will be given a name using this scheme: `featCounts_{species}_{date}.csv`.

## Usage

```
usage: combining_featCount_tables.py [-h] [-d PATH] [-o CUSTOM_FILENAME]

Combines the featureCounts output tables in the target directory.

optional arguments:
  -h, --help            show this help message and exit
  -d PATH, --directory PATH
                        path to target directory. Default: current directory
  -o CUSTOM_FILENAME, --output CUSTOM_FILENAME
                        output filename. Default:
                        featCounts_{species}_{date}.csv
```

## See also

- [fastq-dump_to_featureCounts.sh](fastq-dump_to_featureCounts.md)
