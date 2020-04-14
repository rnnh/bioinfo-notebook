---
layout: default
title: Introduction to the command line
parent: 1. General guides
nav_order: 1
---

# Introduction to the command line

This is a brief introduction to the command line in Ubuntu.
If you are not using Ubuntu, please see the pages with instructions on installing Ubuntu:

- [Windows Subsystem for Linux (WSL, Windows only)](wsl.md)
- [using Ubuntu through a Virtual Machine (Mac or Windows)](ubuntu_virtualbox.md)

## Contents

- [Opening the terminal](#opening-the-terminal)
- [The working directory](#the-working-directory)
- [Cloning the bioinfo-notebook project into your home directory](#cloning-the-bioinfo-notebook-project-into-your-home-directory)
- [Changing working directories](#changing-working-directories)
- [Listing directory content with the `ls` command](#listing-directory-content-with-the-ls-command)
- [Relative paths](#relative-paths)
- [Using the `head` command](#using-the-head-command)
- [The `tail` command](#the-tail-command)
- [Using the `--help` argument](#using-the---help-argument)
- [Running the Linux setup shell script](#running-the-linux-setup-shell-script)
    - [Video demonstration](#video-demonstration)
- [Exercise](#exercise)
- [See also](#see-also)

## Opening the terminal

Once you have installed Ubuntu- either as a virtual machine or a Linux subsystem- open the Ubuntu terminal to use the command line.
If you are using the Ubuntu app from the Microsoft Store in Windows 10 (Windows Subsystem for Linux), you are already using the Ubuntu terminal.
If you are using an Ubuntu virtual machine: click on the "Show Applications" button in Ubuntu (the nine dots in the bottom left corner of the screen), click on the "Type to search..." bar at the top of the screen, type "Terminal" and press Enter to open the command line.

The Ubuntu command line will look like this:

```
(Your UNIX username)@(Your computer's alias):~$ _
```

This is called the *bash prompt*...

- *Your UNIX username* is the username you created when installing Ubuntu
- *Your computer's alias* is the name Ubuntu uses to refer to your computer. This will likely contain the model of your computer (e.g. `Latitude-E7270`).
- The tilde (`~`) indicates that your home directory is the current *working directory*.
- The dollar sign (`$`) indicates that the terminal is using the `bash` shell language.

In examples on this page, `ronan@dell:~ $` will be used as an example bash prompt.

## The working directory

The *working directory* is the directory (also known as "folder") that the command line is currently using.
Any files that you create will be in this directory, and any commands you run will use files in this directory (unless you use a path).

To see the current working directory, type `pwd` into the command line and press Enter (or Return).
This will run the "print working directory" command, which will print the path to the current directory in the terminal.

```bash
ronan@dell:~$ pwd
/home/ronan
```

## Cloning the bioinfo-notebook project into your home directory

We will use the files and directories of the `bioinfo-notebook` project to demonstrate how to use the command line.
This project can be copied into your home directory using the `git clone` command.
The `git clone` command takes the URL of a [GitHub](www.github.com) project, and copies all of the files and directories of that project into the working directory.
To copy the `bioinfo-notebook` project into your Ubuntu system using `git clone`...

1. Copy the URL of this project: <https://github.com/rnnh/bioinfo-notebook>
2. In the Ubuntu command line, type `git clone`. Do *not* press Enter/Return yet.
3. Paste the URL of this project into the command line; either by right-clicking the terminal window and selecting `Paste` in VirtualBox, or just right-clicking the Ubuntu window in WSL.
4. Once `git clone` followed by the URL of this project is in the command line, press Enter/Return.

It should look like this...

```bash
ronan@dell:~$ git clone https://github.com/rnnh/bioinfo-notebook
Cloning into 'bioinfo-notebook'...
remote: Enumerating objects: 178, done.
remote: Counting objects: 100% (178/178), done.
remote: Compressing objects: 100% (138/138), done.
remote: Total 826 (delta 114), reused 68 (delta 40), pack-reused 648
Receiving objects: 100% (826/826), 409.18 KiB | 2.26 MiB/s, done.
Resolving deltas: 100% (505/505), done.
```

The directory `bioinfo-notebook` has now been created in the home directory (`/home/ronan/bioinfo-notebook`).

## Changing working directories

The working directory of the command line can be changed using the `cd` (change directory) command.
Use `pwd` before and after the `cd` command to check which directory you have moved from and to.
Type `cd bioinfo-notebook/` to change the working directory to the `bioinfo-notebook/` directory.

```bash
ronan@dell:~$ pwd
/home/ronan
ronan@dell:~$ cd bioinfo-notebook/
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
```

You do not need to type this command out in full, typing `cd b` and pressing the `Tab` key should complete this command.
Notice that the working directory part of the bash prompt has changed from `~` to `~/bioinfo-notebok`.

When you press `Tab` in the command line, it will try to complete the file or directory name for you.
If there are multiple options, these options will be printed in the terminal.
If there is only one possible option, this option will be filled in.
For example, pressing `Tab` twice after typing `cd D` could give the following directories beginning with "D" as options:

```bash
ronan@dell:~$ cd D
Desktop/   Documents/ Downloads/
```

In Ubuntu (and other Linux systems), directory names are case-sensitive.
This means that `~/downloads/` is a different directory than `~/Downloads/`.

The command `cd ../` can be used to move to move up one directory, `cd ../../` can be used to move up two directories, etc.
In this example, moving up one directory from `/home/ronan/bioinfo-notebook/data/` changes the working directory to `/home/ronan/bioinfo-notebook/`.
Moving up two directories from `/home/ronan/bioinfo-notebook/data/` changes the working directory to `/home/ronan/`.

```bash
ronan@dell:~/bioinfo-notebook/data$ pwd
/home/ronan/bioinfo-notebook/data
ronan@dell:~/bioinfo-notebook/data$ cd ../
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ cd data/
ronan@dell:~/bioinfo-notebook/data$ cd ../../
ronan@dell:~$ pwd
/home/ronan
```
## Listing directory content with the `ls` command

The `ls` command can be used to list the files and directories within the current working directory.
Try using the `ls` command in the `bioinfo-notebook/` directory.

```bash
ronan@dell:~/bioinfo-notebook$ ls
assets  _config.yml  data  docs  envs  LICENSE  README.md  scripts  temp
```

## Relative paths

From the `bioinfo-notebook/` working directory, all of the files in this project can be accessed using *relative paths*.
Without paths, the program only knows which file or directory to look for.
With paths, the program knows which file or directory to look for *and* how to get to it from the working directory.

The `ls` command can be used to list the content of the `bioinfo-notebook/data/` from the `bioinfo-notebook/` working directory using `ls data/`.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ ls data/
example_genome_annotation.gtf  example_nucleotide_sequence.fasta
```

## Using the `head` command

The `head` command can be used to view the first part (the head) of a file or files.
From the `bioinfo-notebook/` working directory, use `head data/example_nucleotide_sequence.fasta` to view the head of that FASTA file.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ head data/example_nucleotide_sequence.fasta 
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACA
CTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTC
CACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTC
ACCCACCGTTACCCTCCAATTACCCATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCA
CCATGACCTACTCACCATACTGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACA
CACACGTGCTTACCCTACCACTTTATACCACCACCACATGCCATACTCACCCTCACTTGTATACTGATTT
TACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTCCACTTCACTC
CATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGG
TCTATACCCTGTGCCATTTACCCATAACGCCCATCATTATCCACATTTTGATATCTATATCTCATTCGGC
```

The `head` command takes optional *arguments*.
Command line arguments are extra instructions given to a program when it runs.
One of the arguments that can be given to `head` is `-n`, which specifies how many lines we want to print per file.
The command `head -n 5 data/example_nucleotide_sequence.fasta` prints the first 5 lines of the file `example_nucleotide_sequence.fasta`.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ head -n 5 data/example_nucleotide_sequence.fasta
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACA
CTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTC
CACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTC
ACCCACCGTTACCCTCCAATTACCCATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCA
```

In the Ubuntu terminal, an asterisk (`*`) acts as a *wildcard*.
This means that any files or directories that can replace this character will replace it.
For example, the `bioinfo-notebook/data/` directory contains two files: `example_genome_annotation.gtf` and `example_nucleotide_sequence.fasta`.
Using the command `head -n 5 data/*` will print the first 5 lines of both of these files.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ head -n 5 data/*
==> data/example_genome_annotation.gtf <==
#gtf-version 2.2
#!genome-build R64
#!genome-build-accession NCBI_Assembly:GCF_000146045.2
#!annotation-source SGD R64-2-1
NC_001133.9	RefSeq	gene	1807	2169	.	-	.	gene_id "YAL068C"; db_xref "GeneID:851229"; gbkey "Gene"; gene "PAU8"; gene_biotype "protein_coding"; locus_tag "YAL068C"; partial "true"; 

==> data/example_nucleotide_sequence.fasta <==
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACA
CTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTC
CACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTC
ACCCACCGTTACCCTCCAATTACCCATATCCAACCCACTGCCACTTACCCTACCATTACCCTACCATCCA
```

## The `tail` command

The `tail` command is the equivalent of the command `head`, but prints the last part of files instead.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ tail -n 5 data/*
==> data/example_genome_annotation.gtf <==
NC_001133.9	RefSeq	exon	225460	226863	.	+	.	gene_id "YAR071W"; transcript_id "NM_001178239.1"; db_xref "GeneID:851299"; gbkey "mRNA"; gene "PHO11"; locus_tag "YAR071W"; partial "true"; product "acid phosphatase PHO11"; exon_number "1"; 
NC_001133.9	RefSeq	CDS	225460	226860	.	+	0	gene_id "YAR071W"; transcript_id "NM_001178239.1"; db_xref "SGD:S000000094"; db_xref "GeneID:851299"; experiment "EXISTENCE:direct assay:GO:0003993 acid phosphatase activity [PMID:8817921]"; gbkey "CDS"; gene "PHO11"; locus_tag "YAR071W"; note "One of three repressible acid phosphatases; glycoprotein that is transported to the cell surface by the secretory pathway; induced by phosphate starvation and coordinately regulated by PHO4 and PHO2; PHO11 has a paralog, PHO12, that arose from a segmental duplication"; product "acid phosphatase PHO11"; protein_id "NP_009434.1"; exon_number "1"; 
NC_001133.9	RefSeq	start_codon	225460	225462	.	+	0	gene_id "YAR071W"; transcript_id "NM_001178239.1"; db_xref "SGD:S000000094"; db_xref "GeneID:851299"; experiment "EXISTENCE:direct assay:GO:0003993 acid phosphatase activity [PMID:8817921]"; gbkey "CDS"; gene "PHO11"; locus_tag "YAR071W"; note "One of three repressible acid phosphatases; glycoprotein that is transported to the cell surface by the secretory pathway; induced by phosphate starvation and coordinately regulated by PHO4 and PHO2; PHO11 has a paralog, PHO12, that arose from a segmental duplication"; product "acid phosphatase PHO11"; protein_id "NP_009434.1"; exon_number "1"; 
NC_001133.9	RefSeq	stop_codon	226861	226863	.	+	0	gene_id "YAR071W"; transcript_id "NM_001178239.1"; db_xref "SGD:S000000094"; db_xref "GeneID:851299"; experiment "EXISTENCE:direct assay:GO:0003993 acid phosphatase activity [PMID:8817921]"; gbkey "CDS"; gene "PHO11"; locus_tag "YAR071W"; note "One of three repressible acid phosphatases; glycoprotein that is transported to the cell surface by the secretory pathway; induced by phosphate starvation and coordinately regulated by PHO4 and PHO2; PHO11 has a paralog, PHO12, that arose from a segmental duplication"; product "acid phosphatase PHO11"; protein_id "NP_009434.1"; exon_number "1"; 
###

==> data/example_nucleotide_sequence.fasta <==
GTGTGGTGATGGATAGTGAGTGGATAGTGAGTGGATGGATGGTGGAGTGGGGGAATGAGACAGGGCATGG
GGTGGTGAGGTAAGTGCCGTGGATTGTGATGATGGAGAGGGAGGGTAGTTGACATGGAGTTAGAATTGGG
TCAGTGTTAGTGTTAGTGTTAGTATTAGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGGTGTG
GGTGTGGGTGTGGGTGTGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGG
```

## Using the `--help` argument

For a lot of command line programs, using the command with the `--help` (or `-h`) argument will display a short help message on how that command is used.

```bash
ronan@dell:~$ head --help
Usage: head [OPTION]... [FILE]...
Print the first 10 lines of each FILE to standard output.
With more than one FILE, precede each with a header giving the file name.

With no FILE, or when FILE is -, read standard input.

Mandatory arguments to long options are mandatory for short options too.
  -c, --bytes=[-]NUM       print the first NUM bytes of each file;
                             with the leading '-', print all but the last
                             NUM bytes of each file
  -n, --lines=[-]NUM       print the first NUM lines instead of the first 10;
                             with the leading '-', print all but the last
                             NUM lines of each file
  -q, --quiet, --silent    never print headers giving file names
  -v, --verbose            always print headers giving file names
  -z, --zero-terminated    line delimiter is NUL, not newline
      --help     display this help and exit
      --version  output version information and exit

NUM may have a multiplier suffix:
b 512, kB 1000, K 1024, MB 1000*1000, M 1024*1024,
GB 1000*1000*1000, G 1024*1024*1024, and so on for T, P, E, Z, Y.

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/head>
or available locally via: info '(coreutils) head invocation'
```

For more in-depth help with a command, the `man` (manual) command is useful, if it is available.
This will open a manual for the program within the terminal.
The `Up` and `Down` arrow keys (or `j` and `k`) can be used to scroll through these manuals, and `q` is used to exit.

```bash
ronan@dell:~$ man head

HEAD(1)                          User Commands                         HEAD(1)

NAME
       head - output the first part of files

SYNOPSIS
       head [OPTION]... [FILE]...

DESCRIPTION
...
 Manual page head(1) line 1 (press h for help or q to quit)
```

## Running the Linux setup shell script

A *shell script* can be used to run commands sequentially, without having to input them individually into the command line.
One of the shell scripts (ending in `.sh`) included in this project is [scripts/linux_setup.sh](linux_setup.md).
This script downloads and installs [conda](conda) and the `bioinfo-notebook` conda environment.
This is a quick way to install command line programs discussed in this project, e.g. [bowtie2](bowtie2.md), [featureCounts](featureCounts.md) and [SAMtools](samtools.md).

A `bash` shell script can be run using `bash` followed by a relative path to the script.
The shell script `linux_setup.sh` also includes help text, which can be accessed with the `--help` or `-h` argument.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ bash scripts/linux_setup.sh --help
linux_setup.sh 
 
 This script downloads and installs Miniconda3, and uses conda to install 
 the 'bioinfo-notebook' virtual environment. 
 
 Before running this script... 
 
 	 1. Please run the following command: 
 	 	 $ sudo apt-get update 
 	 This will ensure that the software installed will be up-to-date. 
 
 	 2. Please ensure that the 'bioinfo-notebook/' directory is in your 
 	 home directory (~). The path to this directory should look like this: 
 	 	 /home/ronan/bioinfo-notebook 
 
 The 'bash' command is used to run this script: 
 	 $ bash scripts/linux_setup.sh 
 
 Optional arguments: 
 	 -h | --help	 show this help text and exit
```

### Video demonstration

In this demonstration, the bioinfo-notebook GitHub project (also known as a repository or repo) is cloned into the home directory of the Linux system (Ubuntu).
This means that all the files for this project will be downloaded from GitHub into the `~/bioinfo-notebook/` directory.
A GitHub repo can be cloned using the command `$ git clone` followed by the URL of the target repo (which can be found on GitHub using the “Clone or download” button).
The Linux setup script is then run from this cloned GitHub repo.

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)

## Exercise

Look at the [structure of the bioinfo-notebook repository](../README.md#repository-structure).
This outlines how this repository (another term for a GitHub project folder) is structured: it outlines which files and directories are in this project.
Most of the files in this project are within subdirectories of the `bioinfo-notebook/` directory.

Once you have read this page, and [copied this project to your Ubuntu system](#cloning-the-bioinfo-notebook-project-into-your-home-directory), try the following small tasks.
These tasks only require one command each.

1. Change the working directory from `bioinfo-notebook/` to `bioinfo-notebook/data/`.
2. Change the working directory from `bioinfo-notebook/data` to `bioinfo-notebook/docs`, using `../` in your command.
3. List the files in the `bioinfo-notebook/docs/` directory.
4. Select a file in the `bioinfo-notebook/docs/` directory, and display the first 6 lines of it using the `head` command.
5. Display the last 2 lines of all the files in the `bioinfo-notebook/docs/` directory, using the `tail` command.
6. From the `bioinfo-notebook/docs/` directory, list the files in the `bioinfo-notebook/envs/` directory.

[Solutions to this exercise](cl_solutions.md).

## See also

- [Windows Subsystem for Linux](wsl.md)
- [Using Ubuntu through a Virtual Machine](ubuntu_virtualbox.md)
- [conda](conda.md)
- [Linux setup script](linux_setup.md)
- [File formats used in bioinformatics](file_formats.md)
- [The DataCamp "Introduction to Shell" interactive course](https://www.datacamp.com/courses/introduction-to-shell-for-data-science)
