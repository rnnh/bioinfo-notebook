---
layout: default
title: OrthoFinder
parent: 2. Program guides
---

# OrthoFinder

OrthoFinder is [a program for phylogenetic orthology inference](https://davidemms.github.io/).
It can be installed using the [orthofinder.yml](../envs/orthofinder.yml) virtual environment using [conda](conda.md).

## Running `OrthoFinder` to find orthologs between sets of FASTA amino acid sequences

`OrthoFinder` can be used to find orthologs between sets of FASTA amino acid files as follows:

```bash
$ orthofinder -t n -S diamond -f path/to/fasta/files/
```

In this command...

1. **`-t`** sets the number of threads/processors to use (*n*).
2. **`-S`** is used to select the search tool OrthoFinder uses. Setting it to [`diamond` is far faster than the default BLAST method](https://github.com/davidemms/OrthoFinder/releases/tag/v2.2.7).
3. **`-f`** is used to select the directory of [FASTA amino acid sequences](file_formats.md#fasta) files you want to compare.

OrthoFinder will create a `Results` directory (ending with the current month and day, e.g. `Results_Sep16/`) in the target directory specified with **`-f`**.
This directory will contain summary statistics of orthologs found between the FASTA files, as well as putative gene duplication events, and phylogenetic trees of the detected orthogroups.

## See also

- [conda](conda.md)
- [File formats used in bioinformatics](file_formats.md)

## Further reading

- [OrthoFinder tutorials](https://davidemms.github.io/menu/tutorials.html)
