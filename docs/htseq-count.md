# htseq-count

Given a file with aligned sequencing reads and a list of genomic features, `htseq-count` can be used to count how many reads map to each feature.

## Aligining reads to a genome annotation using `htseq-count`

`htseq-count` can be used to align reads to a genome annotation as follows:

```
$ htseq-count --format bam sorted_alignment_file.bam genome_annotation > output_file.txt
```

In this command...

1. **`--format`** or **`-f`** is the format of the input data. Possible values are `sam` (for text SAM files) and `bam` (for binary BAM files). Default is `sam`.
A `bam` file is used in this example.
2. **`sorted_alignment_file.bam`** is a `bam` format alignment file, sorted by name.
3. **`genome_annotation`** is the genome annotation file the reads in the `alignment_file` are aligned to (`.gtf` or `.gff`).
4. **`> output_file.txt`** redirects the output (`STDOUT`) to `output_file.txt`.

## The `htseq-count` output file

The program outputs a table with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons.
The names of the special counters all start with a double underscore, to facilitate filtering (**Note:** The double underscore was absent up to version 0.5.4).
The special counters are:

1. **`__no_feature`**: reads (or read pairs) which could not be assigned to any feature (set S as described above was empty).
2. **`__ambiguous`**: reads (or read pairs) which could have been assigned to more than one feature and hence were not counted for any of these, unless the --nonunique all option was used (set S had more than one element).
3. **`__too_low_aQual`**: reads (or read pairs) which were skipped due to the optional minimal alignment quality flag.
4. **`__not_aligned`**: reads (or read pairs) in the SAM/BAM file without an alignment.
5. **`__alignment_not_unique`**: reads (or read pairs) with more than one reported alignment.

## Further reading

1. The `htseq-count` manual: <https://htseq.readthedocs.io/en/release_0.11.1/count.html>
