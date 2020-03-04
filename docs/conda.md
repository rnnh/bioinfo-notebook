# conda

From the website, `conda` provides ["Package, dependency and environment management for any language—Python, R, Ruby, Lua, Scala, Java, JavaScript, C/ C++, FORTRAN, and more"](https://docs.conda.io/en/latest/).

Conda is a package manager allows specific versions of programs to be installed, alongside their dependencies.
Different sets of programs can be installed to different [virtual environments](https://www.anaconda.com/moving-conda-environments/).
A virtual environment is basically a set of programs.
These environments can then be shared, either as a `.yml` file or a `.txt` file.
A `.yml` copy of a conda environment can be used to recreate that environment on another machine, regardless of the operating system platform used.
A `.txt` copy of a conda environment is more explicit: it can be used to create anidentical copy of a conda environment using the same operating system platform as the original machine.
A conda virtual environment is used throughout this project: a [`.yml` copy](../envs/bioinfo-notebook.yml) and an [explicit `.txt` copy](../envs/bioinfo-notebook.txt) of this conda environment are provided.

A conda environment can be activated using the following command...

```
conda activate name_of_environment
```

Once activated, the programs installed in this environment are available.
Conda can be deactivated using `conda deactivate`.

## Demonstration

In this video demonstration, a conda virtual environment is created using [bioinfo-notebook.txt](../envs/bioinfo-notebook.txt).
This virtual environment is then activated using `conda activate bioinfo-notebook`.
Note that the name of the active conda environment is displayed in brackets at the start of the bash prompt: `(name of active environment) ... $`.

[![asciicast](https://asciinema.org/a/305992.svg)](https://asciinema.org/a/305992?autoplay=1)

## Further reading

1. Conda packages: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/packages.html>
2. Conda environments: <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>
