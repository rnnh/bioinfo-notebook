---
layout: default
title: FeatureCounts
parent: 2. Program guides
---


# FeatureCounts

`featureCounts` is a program that counts how many reads map to genomic features, such as genes, exon, promoter and genomic bins.

## Counting how many reads align to each gene in a genome annotation using `featureCounts`

`featureCounts` can be used to count how many reads align to genes as follows:

```
$ featureCounts -p -O -T n -a example_genome_annotation.gtf -o example_featureCounts_output.txt sorted_example_alignment.bam
```

In this command...

1. **`-p`** species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
2. **`-O`** assigns reads to all their overlapping meta-features.
3. **`-T`** specifies the number (*`n`*) of threads to be used.
4. **`-a`** is the genome annotation file (`example_genome_annotation.gtf`).
5. **`-o`** specifies the name of the output file, which includes the read counts (`example_featureCounts_output.txt`).
6. **`sorted_example_alignment.bam`** is an alignment file: in this file, the reads we want to count are aligned to the same genome as the annotation file.

### Demonstration

In this video, `featureCounts` is used to assign reads in an alignment file (`sorted_example_alignment.bam`) to genes in a genome annotation file (`example_genome_annotation.gtf`).

[![asciicast](https://asciinema.org/a/306584.svg)](https://asciinema.org/a/306584?autoplay=1)

## More important options for `featureCounts`

1. **`-s`** specifies strand-specific read counting. `0` for unstranded reads, `1` for stranded reads and `2` for reversely stranded reads. This depends on the library used in the sequencing protocol.

## Further reading

1. The `subread` user guide: <http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf>
2. The `featureCounts` paper: <https://doi.org/10.1093/bioinformatics/btt656>
