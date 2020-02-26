# fastq-dump

```
fastq-dump --gzip --skip-technical  --readids --read-filter pass 
--dumpbase --split-3 --clip --outdir path/to/reads/ SRR_ID
```

- **`-O | --outdir <path>`** *(Optional)* Output directory, default is current working directory
 ('.').
- **`--gzip`** Compress output using gzip. Gzip archived reads can be read directly by bowtie2.
- **`--skip-technical`** Dump only biological reads.
- **`-I | --readids`** Append read id after spot id as 'accession.spot.readid' on 
defline. With this flag, one sequence gets appended the ID .1 and the other .2.
- **`--read-filter pass`** Only returns reads that pass filtering (without `N`s).
- **`-B | --dumpbase`** Formats sequence using base space (default for other than 
SOLiD). Included to avoid colourspace (in which pairs of bases are represented by numbers).
- **`-W|--clip`** Some of the sequences in the SRA contain tags that are used e.g. for whole genome amplification and need to be removed. This will remove those sequence.
- **`--split-3`** separates the read into left and right ends, but if there is a left end without a matching right end, or a right end without a matching left end, they will be put in a single file.
- **`--split-files`** splits the FASTQ reads into two files: one file for mate 1s (`...1`), and another for mate 2s (`..._2`).
- **`--split-spot`** splits the FASTQ reads into two (mate 1s and mate 2s) within one file. "split-spot will give you an 8-line fastq format where forward precedes reverse" (<https://www.biostars.org/p/178586/>).
- **Reference:** <https://edwards.sdsu.edu/research/fastq-dump/>
