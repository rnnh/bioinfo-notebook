# htseq-count

<https://htseq.readthedocs.io/en/release_0.11.1/count.html#count>

`htseq-count [options] <alignment_files> <gff_file>`

The program outputs a table with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons. The names of the special counters all start with a double underscore, to facilitate filtering. (Note: The double underscore was absent up to version 0.5.4). The special counters are:

- `__no_feature`: reads (or read pairs) which could not be assigned to any feature (set S as described above was empty).
- `__ambiguous`: reads (or read pairs) which could have been assigned to more than one feature and hence were not counted for any of these, unless the --nonunique all option was used (set S had more than one element).
- `__too_low_aQual`: reads (or read pairs) which were skipped due to the -a option, see below
- `__not_aligned`: reads (or read pairs) in the SAM file without alignment
- `__alignment_not_unique`: reads (or read pairs) with more than one reported alignment. These reads are recognized from the NH optional SAM field tag. (If the aligner does not set this field, multiply aligned reads will be counted multiple times, unless they getv filtered out by due to the -a option.) Note that if the --nonunique all option was used, these reads (or read pairs) are still assigned to features.

## Options

- `-f <format>, --format=<format>`

Format of the input data. Possible values are `sam` (for text SAM files) and `bam` (for binary BAM files). Default is `sam`.

- `r <order>, --order=<order>`

For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the `samtools` sort function of `samtools` to sort it. Use this option, with `name` or `pos` for `<order>` to indicate how the input data has been sorted. The default is `name`.

If `name` is indicated, `htseq-count` expects all the alignments for the reads of a given read pair to appear in adjacent records in the input data. For `pos`, this is not expected; rather, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found. While, strictly speaking, the latter will also work with unsorted data, sorting ensures that most alignment mates appear close to each other in the data and hence the buffer is much less likely to overflow.

- `--max-reads-in-buffer=<number>`

When `<alignment_file>` is paired end sorted by position, allow only so many reads to stay in memory until the mates are found (raising this number will use more memory). Has no effect for single end or paired end sorted by name. (default: 30000000)

- `-s <yes/no/reverse>, --stranded=<yes/no/reverse>`

Whether the data is from a strand-specific assay (default: `yes`).

For `stranded=no`, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For `stranded=yes` and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For `stranded=reverse`, these rules are reversed.

- `-a <minaqual>, --a=<minaqual>`

skip all reads with alignment quality lower than the given minimum value (default: 10 — Note: the default used to be 0 until version 0.5.4.)

- `-t <feature type>, --type=<feature type>`

feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon)

- `-i <id attribute>, --idattr=<id attribute>`

GFF attribute to be used as feature ID. Several GFF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is `gene_id`.

- `--additional-attr=<id attributes>`

Additional feature attributes, which will be printed as an additional column after the primary attribute column but before the counts column(s). The default is none, a suitable value to get gene names using an Ensembl GTF file is `gene_name`.
To use more than one additional attribute, repeat the option in the command line more than once, with a single attribute each time, e.g. `--additional-attr=gene_name` `--additional_attr=exon_number`.

- `-m <mode>, --mode=<mode>`

Mode to handle reads overlapping more than one feature. Possible values for `<mode>` are `union`, `intersection-strict` and `intersection-nonempty` (default: `union`)

- `--nonunique=<nonunique mode>`

Mode to handle reads that align to or are assigned to more than one feature in the overlap `<mode>` of choice (see `-m` option). `<nonunique mode>` are none and all (default: none)

- `--secondary-alignments=<mode>`

Mode to handle secondary alignments (SAM flag 0x100). <mode> can be score and ignore (default: score)

- `--supplementary-alignments=<mode>`

Mode to handle supplementary/chimeric alignments (SAM flag 0x800). `<mode>` can be score and ignore (default: score)

- `-o <samout>, --samout=<samout>`

write out all SAM alignment records into an output SAM file called `<samout>`, annotating each line with its assignment to a feature or a special counter (as an optional field with tag ‘XF’)

- `-q, --quiet`

suppress progress report and warnings

- `-h, --help`

Show a usage summary and exit
