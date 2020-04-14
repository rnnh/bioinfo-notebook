---
layout: default
title: Command line exercise solutions
nav_exclude: true
---

# Command line exercise solutions

The `pwd` commands in these solutions are added to clarify the working directories used.

1. Change the working directory from `bioinfo-notebook/` to `bioinfo-notebook/data/`.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ cd data/
ronan@dell:~/bioinfo-notebook/data$ pwd
/home/ronan/bioinfo-notebook/data
```

2. Change the working directory from `bioinfo-notebook/data` to `bioinfo-notebook/docs`, using `../` in your command.

```bash
ronan@dell:~/bioinfo-notebook/data$ pwd
/home/ronan/bioinfo-notebook/data
ronan@dell:~/bioinfo-notebook/data$ cd ../docs/
ronan@dell:~/bioinfo-notebook/docs$ pwd
/home/ronan/bioinfo-notebook/docs
```

3. List the files in the `bioinfo-notebook/docs/` directory.

```bash
ronan@dell:~/bioinfo-notebook/docs$ pwd
/home/ronan/bioinfo-notebook/docs
ronan@dell:~/bioinfo-notebook/docs$ ls
bowtie2.md                      file_formats.md
bowtie.md                       htseq-count.md
cl_intro.md                     linux_setup.md
cl_solutions.md                 part1.md
combining_featCount_tables.md   part2.md
conda.md                        part3.md
fasterq-dump.md                 samtools.md
fastq-dump.md                   to_do.md
fastq-dump_to_featureCounts.md  ubuntu_virtualbox.md
featureCounts.md                wsl.md
```

4. Select a file in the `bioinfo-notebook/docs/` directory, and display the first 6 lines of it using the `head` command.

```bash
ronan@dell:~/bioinfo-notebook/docs$ pwd
/home/ronan/bioinfo-notebook/docs
ronan@dell:~/bioinfo-notebook/docs$ head cl_solutions.md 
---
layout: default
title: Command line exercise solutions
nav_exclude: true
---

# Command line exercise solutions

1. Change the working directory from `bioinfo-notebook/` to `bioinfo-notebook/data/`.
```

5. Display the last 2 lines of all the files in the `bioinfo-notebook/docs/` directory, using the `tail` command.

```bash
ronan@dell:~/bioinfo-notebook/docs$ pwd
/home/ronan/bioinfo-notebook/docs
ronan@dell:~/bioinfo-notebook/docs$ tail -n 2 *
==> bowtie2.md <==

1. The `bowtie2` manual: <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>

==> bowtie.md <==

1. The `bowtie` manual: <http://bowtie-bio.sourceforge.net/manual.shtml>

==> cl_intro.md <==
- [File formats used in bioinformatics](file_formats.md)
- [The DataCamp "Introduction to Shell" interactive course](https://www.datacamp.com/courses/introduction-to-shell-for-data-science)

==> cl_solutions.md <==
5. Display the last 2 lines of all the files in the `bioinfo-notebook/docs/` directory, using the `tail` command.
6. From the `bioinfo-notebook/docs/` directory, list the files in the `bioinfo-notebook/envs/` directory.

==> combining_featCount_tables.md <==

- [fastq-dump_to_featureCounts.sh](fastq-dump_to_featureCounts.md)

==> conda.md <==
2. Conda packages: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/packages.html>
3. Conda environments: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>

==> fasterq-dump.md <==

1. [How to use fasterq-dump from the sra-tools wiki on GitHub](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)

==> fastq-dump.md <==

1. Rob Edward's notes on `fastq-dump`: <https://edwards.sdsu.edu/research/fastq-dump/>

==> fastq-dump_to_featureCounts.md <==

1. [fastq-dump_to_featureCounts.sh on GitHub](https://github.com/rnnh/bioinfo-notebook/blob/master/scripts/fastq-dump_to_featureCounts.sh)

==> featureCounts.md <==
1. The `subread` user guide: <http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf>
2. The `featureCounts` paper: <https://doi.org/10.1093/bioinformatics/btt656>

==> file_formats.md <==
- [GTF2.2: A Gene Annotation Format (Revised Ensembl GTF)](http://mblab.wustl.edu/GTF22.html)
- [GFF3 Specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

==> htseq-count.md <==

1. The `htseq-count` manual: <https://htseq.readthedocs.io/en/release_0.11.1/count.html>

==> linux_setup.md <==
- [Using Ubuntu through a Virtual Machine](ubuntu_virtualbox.md)
- [Windows Subsystem for Linux](wsl.md)

==> part1.md <==

These are general guides for installing Ubuntu, using the command line, and the types of files used in bioinformatics.

==> part2.md <==

These are guides to individual programs.

==> part3.md <==

These are scripts that use the programs discussed in this project.

==> samtools.md <==
- [Alignment formats](file_formats.md#alignment-formats)
- The `samtools` manual: <https://www.htslib.org/doc/samtools.html>

==> to_do.md <==
- Add page on `trimmomatic`
- Entry on BED/bigWig

==> ubuntu_virtualbox.md <==
- [What is a Virtual Machine?](https://azure.microsoft.com/en-us/overview/what-is-a-virtual-machine/)
- [How to Install Ubuntu on VirtualBox](https://www.wikihow.com/Install-Ubuntu-on-VirtualBox)

==> wsl.md <==
- [Using Ubuntu through a Virtual Machine](ubuntu_virtualbox.md) 
- [conda](conda.md)
```

6. From the `bioinfo-notebook/docs/` directory, list the files in the `bioinfo-notebook/envs/` directory.

```bash
ronan@dell:~/bioinfo-notebook/docs$ pwd
/home/ronan/bioinfo-notebook/docs
ronan@dell:~/bioinfo-notebook/docs$ ls ../envs/
bioinfo-notebook.txt  bioinfo-notebook.yml
```
