# Linux setup script

[linux_setup.sh](../scripts/linux_setup.md) is a `bash` shell script that...

1. Downloads and installs [Miniconda3](conda.md)
2. Installs the `bioinfo-notebook` [virtual environment using conda](conda.md#cloning-and-activating-a-conda-environment)

This will use around 2.3 GB of hard disk space in total.

If you are using a Linux system that does not have Anaconda/Miniconda installed, this script will set up everything you need to follow the guides on this website.
If you are using a freshly installed [Ubuntu virtual machine](ubuntu_virtualbox.md) or [Ubuntu through Windows Subsystem for Linux](wsl.md), this script is the ideal way to set up your new system.

## Demonstration

This is a video demonstration of [linux_setup.sh](../scripts/linux_setup.sh).

In this demonstration, the [bioinfo-notebook GitHub repository](https://github.com/rnnh/bioinfo-notebook) (or "repo") is cloned into the home directory of the Linux system (Ubuntu).
This means that all the files for this project will be downloaded from GitHub into the `~/bioinfo-notebook/` directory.
A GitHub repo can be cloned using the command `$ git clone` followed by the URL of the target repo (which can be found on GitHub using the "Clone or download" button).
The Linux setup script is then run from this cloned GitHub repo.

[![asciicast](https://asciinema.org/a/314853.svg)](https://asciinema.org/a/314853?autoplay=1)

## Usage

```
This script downloads and installs Miniconda3, and uses conda to install \n
the 'bioinfo-notebook' virtual environment.

Before running this script...

	1. Please run the following command:
		$ sudo apt-get update
	This will ensure that the software installed will be up-to-date.

	2. Please ensure that the 'bioinfo-notebook/' directory is in your
	home directory (~). The path to this directory should look like this:
		$HOME/bioinfo-notebook

The 'bash' command is used to run this script:
	$ bash $0

Optional arguments:
	-h | --help	show this help text and exit
```
