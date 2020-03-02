# samtools

`samtools` is a set of utilities that manipulate alignments in the BAM format.
It imports from and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to retrieve reads in any regions swiftly. 

## Converting a `sam` alignment file to a sorted, indexed `bam` file using `samtools`

Sequence Alignment Map (SAM/`.sam`) is a text-based file is a text-based file format for sequence alignments.
It's binary equivalent is Binary Alignment Map (BAM/`.bam`), which stores the same data as a compressed binary file.
A binary file for a sequence alignment is preferable over a text file, as binary files are faster to work with.
A SAM alignment file (`example_alignment.sam`) can be converted to a BAM alignment using `samtools view`.

```
$ samtools view -@ n -Sb -o example_alignment.bam example_alignment.sam
```

In this command...

1. **`-@`** sets the number (*`n`*) of threads/CPUs to be used. This flag is optional and can be used with other `samtools` commands.
2. **`-Sb`** specifies that the input is in SAM format (`S`) and the output will be be BAM format(`b`).
3. **`-o`** sets the name of the output file (`example_alignment.bam`).
4. **`example_alignment.sam`** is the name of the input file.

Now that the example alignment is in BAM format, we can sort it using `samtools sort`.
Sorting this alignment will allow us to create a index.

```
$ samtools sort -O bam -o sorted_example_alignment.bam example_alignment.bam
```

In this command...

1. **`-O`** specifies the output format (`bam`, `sam`, or `cram`).
2. **`-o`** sets the name of the output file (`sorted_example_alignment.bam`).
3. **`example_alignment.bam`** is the name of the input file.

This sorted BAM alignment file can now be indexed using `samtools index`.
Indexing speeds allows fast random access to this alignment, allowing the information in the alignment file to be processed faster.

```
$ samtools index sorted_example_alignment.bam
```

In this command...

1. **`sorted_example_alignment.bam`** is the name of the input file.

## Simulating short reads using `samtools wgsim`

`wgsim` is a `samtools` program that can simulate short sequencing reads from a reference genome.
This is useful for creating FASTQ files to practice with.

```
$ wgsim example_nucleotide_sequence.fasta example_reads_1.fastq example_reads_2.fastq
```

In this command...

1. **`example_nucleotide_sequence.fasta`** is the reference genome input.
2. **`example_reads_1.fastq`** and **`example_reads_2.fastq`** are the names of the simulated read output files.

### Demonstration

In this video, `wgsim` is used to simulate reads from `example_nucleotide_sequence.fasta`.


<script id="asciicast-m89gXtx4cKRnKpI6amWj3BEAH" src="https://asciinema.org/a/m89gXtx4cKRnKpI6amWj3BEAH.js" async></script>

## Further reading

1. The `samtools` manual: <https://www.htslib.org/doc/samtools.html>
