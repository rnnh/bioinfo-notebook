# `samtools`

<https://www.htslib.org/doc/samtools.html>

## `view`

`samtools view [options] in.sam|in.bam|in.cram [region...]`

With no options or regions specified, prints all alignments in the specified input alignment file (in SAM, BAM, or CRAM format) to standard output in SAM format (with no header). 

The `-b`, `-C`, `-1`, `-u`, `-h`, `-H`, and `-c` options change the output format from the default of headerless SAM, and the `-o` and `-U` options set the output file name(s).

The `-t` and `-T` options provide additional reference data. One of these two options is required when SAM input does not contain @SQ headers, and the `-T` option is required whenever writing CRAM output.

The `-L`, `-M`, `-r`, `-R`, `-s`, `-q`, `-l`, `-m`, `-f`, `-F`, and `-G` options filter the alignments that will be included in the output to only those alignments that match certain criteria.

The `-x` and `-B` options modify the data which is contained in each alignment.

Finally, the `-@` option can be used to allocate additional threads to be used for compression, and the `-?` option requests a long help message. 

- `-b` Output in the BAM format. 
- `-o` *FILE* Output to *FILE [stdout]*.
- `-@` *INT* Number of BAM compression threads to use in addition to main thread [0]. 
- `-S` Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input. 

## `sort`

`samtools sort [-l level] [-m maxMem] [-o out.bam] [-O format] [-n] [-t tag]`
`[-T tmpprefix] [-@ threads] [in.sam|in.bam|in.cram]`

Sort alignments by leftmost coordinates, or by read name when `-n` is used. An appropriate `@HD-SO` sort order header tag will be added or an existing one updated if necessary.


## `index`

`samtools index [-bc] [-m INT] aln.bam|aln.cram [out.index]`

Index a coordinate-sorted BAM or CRAM file for fast random access. (Note that this does not work with SAM files even if they are bgzip compressed — to index such files, use tabix(1) instead.)

This index is needed when region arguments are used to limit samtools view and similar commands to particular regions of interest.

If an output filename is given, the index file will be written to out.index. Otherwise, for a CRAM file aln.cram, index file aln.cram.crai will be created; for a BAM file aln.bam, either aln.bam.bai or aln.bam.csi will be created, depending on the index format selected.
