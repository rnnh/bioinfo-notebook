---
layout: default
title: Introduction to the command line
parent: 1. General guides
nav_order: 1
---

# Introduction to the command line

This is an introduction to the *UNIX (Mac/Linux)* command line for absolute beginners.
It covers some commands and concepts that are widely used in the UNIX command line.
Examples are provided throughout.
It is recommend that you try out commands as you read about them in your own UNIX system, as well as the exercise towards the end of the page.

If you are using Windows, please see these pages with instructions on installing Ubuntu (a UNIX operating system):

- [Windows Subsystem for Linux (WSL, Windows only)](wsl.md)
- [Using Ubuntu through a Virtual Machine (Mac or Windows)](ubuntu_virtualbox.md)

## Contents

- [Opening the terminal](#opening-the-terminal)
- [The working directory](#the-working-directory)
- [Cloning the `bioinfo-notebook` project into your home directory](#cloning-the-bioinfo-notebook-project-into-your-home-directory)
- [Changing working directories](#changing-working-directories)
- [Comments and broken lines](#comments-and-broken-lines)
- [Listing directory content with the `ls` command](#listing-directory-content-with-the-ls-command)
- [Relative paths](#relative-paths)
- [Making directories with `mkdir`](#making-directories-with-mkdir)
- [Removing files and directories with the `rm` command](#removing-files-and-directories-with-the-rm-command)
- [Using the `head` command](#using-the-head-command)
- [The `tail` command](#the-tail-command)
- [Using the `--help` argument](#using-the---help-argument)
- [Downloading files with `wget`](#downloading-files-with-wget)
- [Moving and copying with `mv` and `cp`](#moving-and-copying-with-mv-and-cp)
- [Running the Linux setup shell script](#running-the-linux-setup-shell-script)
    - [Video demonstration](#video-demonstration)
- [Exercise](#exercise)
- [See also](#see-also)

## Opening the terminal

Once you are using a UNIX operating system (i.e. a Mac system, a Linux system, or Ubuntu through either a virtual machine or a Linux subsystem on a Windows machine) open the *terminal* to use the command line.
The terminal is the window in which the command line runs.

- If you are using a Mac, click on Spotlight Search (the magnifying glass icon in the top-right corner of the screen), type "Terminal", and open the Terminal application.
- If you are using the Ubuntu app from the Microsoft Store in Windows 10 (Windows Subsystem for Linux), you are already using the Ubuntu terminal.
- If you are using an Ubuntu virtual machine: click on the "Show Applications" button in Ubuntu (the nine dots in the bottom left corner of the screen), click on the "Type to search..." bar at the top of the screen, type "Terminal" and press `Enter` to open the command line.

The UNIX command line will look like this:

```
(Your UNIX username)@(Your computer's alias):~$ _
```

This is called the *bash prompt*...

- *Your UNIX username* is the your system username.
- *Your computer's alias* is the name UNIX uses to refer to your computer. This will likely contain the model of your computer (e.g. `Latitude-E7270`).
- The tilde (`~`) indicates that your home directory is the current *working directory*. The home directory is located at `/home/` followed by your UNIX username.
- The dollar sign (`$`) indicates that the terminal is using the `bash` shell language.

In examples on this page, `ronan@dell:~ $` will be used as an example bash prompt.

## The working directory

A *directory* is the same as a folder on your desktop.
If you have a "Pictures" folder on your computer's desktop, this folder is a directory within the desktop directory: its [path](#relative-paths) is `Desktop/Pictures/`.
The *working directory* is the directory that the command line is currently using.
Any files that you create will be in this directory, and any commands you run will use files in this directory (unless you use a [path](#relative-paths)).

To see the current working directory, type `pwd` into the command line and press `Enter` (or `Return`).
This will run the "print working directory" command, which will print the path to the current directory in the terminal.
In the context of command line programs, "print" means "display in the terminal".

```bash
ronan@dell:~$ pwd
/home/ronan
```

## Cloning the bioinfo-notebook project into your home directory

We will use the files and directories of the `bioinfo-notebook` project to demonstrate how to use the command line.
This project can be copied into your home directory using the `git clone` command.
The `git clone` command takes the URL of a [GitHub](https://github.com/) project, and copies all of the files and directories of that project into the working directory.
To copy the `bioinfo-notebook` project into your UNIX system using `git clone`...

1. Copy the URL of this project: <https://github.com/rnnh/bioinfo-notebook>
2. In the UNIX command line, type `git clone`. Do *not* press `Enter`/`Return` yet.
3. Paste the URL of this project into the command line; either by right-clicking the terminal window and selecting `Paste` in VirtualBox, or just right-clicking the Ubuntu window in WSL.
4. Once `git clone` followed by the URL of this project is in the command line, press `Enter`/`Return`.

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

In UNIX, directory names are case-sensitive.
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

## Comments and broken lines

The number sign or hash (`#`) is used for *comments* in bash.
Comments are used in programming to explain or annotate code; they are not run as code.
Anything written after a `#` symbol in the command line is a comment.

```bash
ronan@dell:~/bioinfo-notebook$ pwd # This prints the working directory
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ ls # This lists the contents of the directory
assets  _config.yml  data  docs  envs  LICENSE  README.md  scripts  temp
```

Long lines of code can be hard to read.
For legibility, long lines of command line code can be broken up with the backslash (`\`).
The `\` characters split lines to make them more readable, they do not affect how the code functions.

```bash
ronan@dell:~/bioinfo-notebook/data$ head example_nucleotide_sequence.fasta 
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
ronan@dell:~/bioinfo-notebook/data$ head \
> example_nucleotide_sequence.fasta 
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

## Listing directory content with the `ls` command

The `ls` command can be used to list the files and directories within the current working directory.
Try using the `ls` command in the `bioinfo-notebook/` directory.

```bash
ronan@dell:~/bioinfo-notebook$ ls
assets  _config.yml  data  docs  envs  LICENSE  README.md  scripts  temp
```

## Relative paths

From the `bioinfo-notebook/` working directory, all of the files in this project can be accessed using *relative paths*.
Without paths, the program only knows which file or directory to look for in the working directory.
With paths, the program knows which file or directory to look for *and* how to get to it from the working directory.

The `ls` command can be used to list the content of the `bioinfo-notebook/data/` from the `bioinfo-notebook/` working directory using `ls data/`.

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ ls data/
example_genome_annotation.gtf  example_nucleotide_sequence.fasta
```

## Making directories with `mkdir`

The command `mkdir` can be used to create a directory.
Here is the [--help](#using-the---help-argument) text for the command `mkdir`:

```bash
ronan@dell:~$ mkdir --help
Usage: mkdir [OPTION]... DIRECTORY...
Create the DIRECTORY(ies), if they do not already exist.

Mandatory arguments to long options are mandatory for short options too.
  -m, --mode=MODE   set file mode (as in chmod), not a=rwx - umask
  -p, --parents     no error if existing, make parent directories as needed
  -v, --verbose     print a message for each created directory
  -Z                   set SELinux security context of each created directory
                         to the default type
      --context[=CTX]  like -Z, or if CTX is specified then set the SELinux
                         or SMACK security context to CTX
      --help     display this help and exit
      --version  output version information and exit

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/mkdir>
or available locally via: info '(coreutils) mkdir invocation'
```

In this example, `mkdir` is used to create a directory called `example`, which is located in (and is therefore a subdirectory of) the home directory (`~`).
This is the equivalent of making a new folder within a folder on your desktop.

```bash
ronan@dell:~$ clear
ronan@dell:~$ pwd
/home/ronan
ronan@dell:~$ mkdir example
ronan@dell:~$ cd example/
ronan@dell:~/example$ pwd
/home/ronan/example
```

## Removing files and directories with the `rm` command

The `rm` (remove) command is used to remove files and directories.
This command must be used *with care,* as it is not the same as deleting a file or directory on a desktop.
Deleted files on a desktop are moved to a Recycle Bin or Trash folder, and are only permanently deleted once they are removed from this folder.
When using the `rm` command, files and directories are permanently deleted; they cannot be recovered.

Here is the [--help](#using-the---help-argument) text for the `rm` command:

```bash
ronan@dell:~$ rm --help
Usage: rm [OPTION]... [FILE]...
Remove (unlink) the FILE(s).

  -f, --force           ignore nonexistent files and arguments, never prompt
  -i                    prompt before every removal
  -I                    prompt once before removing more than three files, or
                          when removing recursively; less intrusive than -i,
                          while still giving protection against most mistakes
      --interactive[=WHEN]  prompt according to WHEN: never, once (-I), or
                          always (-i); without WHEN, prompt always
      --one-file-system  when removing a hierarchy recursively, skip any
                          directory that is on a file system different from
                          that of the corresponding command line argument
      --no-preserve-root  do not treat '/' specially
      --preserve-root   do not remove '/' (default)
  -r, -R, --recursive   remove directories and their contents recursively
  -d, --dir             remove empty directories
  -v, --verbose         explain what is being done
      --help     display this help and exit
      --version  output version information and exit

By default, rm does not remove directories.  Use the --recursive (-r or -R)
option to remove each listed directory, too, along with all of its contents.

To remove a file whose name starts with a '-', for example '-foo',
use one of these commands:
  rm -- -foo

  rm ./-foo

Note that if you use rm to remove a file, it might be possible to recover
some of its contents, given sufficient expertise and/or time.  For greater
assurance that the contents are truly unrecoverable, consider using shred.

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/rm>
or available locally via: info '(coreutils) rm invocation'
```

To remove files, simply type `rm` followed by the name of the file, and press `Enter`/`Return`:

```bash
ronan@dell:~/example$ pwd
/home/ronan/example
ronan@dell:~/example$ ls
example_file.txt
ronan@dell:~/example$ rm example_file.txt # Remove the file 'example_file.txt'
ronan@dell:~/example$ ls
ronan@dell:~/example$ 
```

To remove directories, `rm -r` must be typed, followed by the name of the directory.
**Use the `rm -r` command with caution,** as it will permanently delete entire directories and *all* of their contents.
This includes files and other directories within this directory.

In this example, `rm -r` is used to remove the `example/` directory from the home directory:

```bash
ronan@dell:~/example$ pwd
/home/ronan/example
ronan@dell:~/example$ cd ../
ronan@dell:~$ pwd
/home/ronan
ronan@dell:~$ rm -r example/
ronan@dell:~$ cd example
bash: cd: example: No such file or directory
ronan@dell:~$ 
```

## Using the `head` command

The `head` command can be used to view the first part (the head) of a file or files.
This command is useful for examining very large files quickly.
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

In the UNIX terminal, an asterisk (`*`) acts as a *wildcard*.
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

The asterisk command is especially useful for selecting files with the same *file extension*.
The file extension is the part of the filename after the full stop that specifies the file type: for example, a file ending in `.txt` is a text file.
In a directory with many text files, `*.txt` selects all of the text files.

## The `tail` command

The `tail` command is the equivalent of the command `head`, but prints the last part of files instead.
This is useful for quickly examining the information at the ends of files.

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

## Downloading files with `wget`

The `wget` command can be used to download files from the directory using the command line.
It can take a URL as an argument, and download the file found at that URL into the working directory.

In this example, `wget` is used to download the amino acid sequence of [haemoglobin subunit alpha from UniProt](https://www.uniprot.org/uniprot/P69905):

```bash
ronan@dell:~/bioinfo-notebook$ pwd
/home/ronan/bioinfo-notebook
ronan@dell:~/bioinfo-notebook$ wget \
> https://www.uniprot.org/uniprot/P69905.fasta
--2020-04-15 16:52:27--  https://www.uniprot.org/uniprot/P69905.fasta
Resolving www.uniprot.org (www.uniprot.org)... 128.175.245.202, 193.62.192.81
Connecting to www.uniprot.org (www.uniprot.org)|128.175.245.202|:443... connected.
HTTP request sent, awaiting response... 200 
Length: 233 [text/plain]
Saving to: ‘P69905.fasta’

P69905.fasta        100%[===================>]     233  --.-KB/s    in 0s      

2020-04-15 16:52:28 (7.91 MB/s) - ‘P69905.fasta’ saved [233/233]

ronan@dell:~/bioinfo-notebook$ head P69905.fasta # Examining the start of the downloaded file 
>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
```

## Moving and copying with `mv` and `cp`

The `mv` and `cp` commands can be used to move and copy files or directories from the command line.
Both of these commands take a target file/directory and a destination file/directory as arguments.
The `cp` commands creates a copy of a file/directory, whereas `mv` moves it to a new destination.
The `mv` command can be used to rename files, by moving it to the same directory, but with a new file name.

```bash
ronan@dell:~/bioinfo-notebook$ mv P69905.fasta data/ # Moving the P69905.fasta file to the data/ directory

ronan@dell:~/bioinfo-notebook$ cd data/ # Changing directory to data/

ronan@dell:~/bioinfo-notebook/data$ ls # Listing files in data/ to confirm that the FASTA file has been moved
example_genome_annotation.gtf      P69905.fasta
example_nucleotide_sequence.fasta  S_cere_GCF_000146045.2_R64_genomic.fna

ronan@dell:~/bioinfo-notebook/data$ mv P69905.fasta haem.fasta # Using mv to rename the FASTA file to haem.fasta

ronan@dell:~/bioinfo-notebook/data$ ls # Listing the files in the directory
example_genome_annotation.gtf      haem.fasta
example_nucleotide_sequence.fasta  S_cere_GCF_000146045.2_R64_genomic.fna

ronan@dell:~/bioinfo-notebook/data$ cp haem.fasta ../ # Creating a copy of haem.fasta in the directory above data/, i.e. the bioinfo-notebook/ directory

ronan@dell:~/bioinfo-notebook/data$ ls # Listing the files in the directory
example_genome_annotation.gtf      haem.fasta
example_nucleotide_sequence.fasta  S_cere_GCF_000146045.2_R64_genomic.fna

ronan@dell:~/bioinfo-notebook/data$ cd ../ # Changing directory to bioinfo-notebook/

ronan@dell:~/bioinfo-notebook$ ls # Listing the files in the directory
assets       data  envs        LICENSE    scripts
_config.yml  docs  haem.fasta  README.md  temp
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

In this demonstration, the bioinfo-notebook GitHub project (also known as a repository or repo) is cloned into the home directory of the UNIX system (in this case, the UNIX system used is Ubuntu).
This means that all the files for this project will be downloaded from GitHub into the `~/bioinfo-notebook/` directory.
A GitHub repo can be cloned using the command `$ git clone` followed by the URL of the target repo (which can be found on GitHub using the “Clone or download” button).
The Linux setup script is then run from this cloned GitHub repo.

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)

## Exercise

Look at the [structure of the bioinfo-notebook repository](../README.md#repository-structure).
This outlines how this repository (another term for a GitHub project folder) is structured: it outlines which files and directories are in this project.
Most of the files in this project are within subdirectories of the `bioinfo-notebook/` directory.

Once you have read this page, and [copied this project to your UNIX system](#cloning-the-bioinfo-notebook-project-into-your-home-directory), try the following small tasks.
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
