# fastq-dump to featureCounts.sh

[fastq-dump_to_featureCounts.sh](../scripts/fastq-dump_to_featureCounts.sh) is a `bash` script that...

1. Downloads FASTQ reads from NCBI's SRA using [fastq-dump](fastq-dump.md)
2. Indexes a reference genome and aligns reads to that index using [bowtie2](bowtie2.md)
3. Converts the alignment file created by bowtie2 to BAM format and sorts it using [samtools](samtools.md)
4. Assigns the read alignments to genes in a genome annotation file using [featureCounts](featureCounts.md)

## Demonstration

This is a video demonstration of [fastq-dump_to_featureCounts.sh](../scripts/fastq-dump_to_featureCounts.sh).

During this demonstration, the full genome sequence and genome annotation for [*Saccharomyces cerevisiae* S288C](https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2) are used. The files [example_nucleotide_sequence.fasta](../data/example_nucleotide_sequence.fasta) and [example_genome_annotation.gtf](../data/example_genome_annotation.gtf) are fragments of the nucleotide sequence and annotation for this genome. [RNA-Seq reads for *Saccharomyces cerevisiae* (SRR8933535)](https://www.ncbi.nlm.nih.gov/sra/SRR8933535) are used as the example FASTQ files in this demonstration.

[![asciicast](https://asciinema.org/a/307425.svg)](https://asciinema.org/a/307425?autoplay=1)

## See also

1. [fastq-dump_to_featureCounts.sh on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/fastq-dump_to_featureCounts.sh)
