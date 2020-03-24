#! /bin/bash/

'''
Linux Setup Script

This script downloads and installs Miniconda3, and uses conda to install
the "bioinfo-notebook" virtual environment
'''

echo Checking if bioinfo-notebook/ is in the home directory...
sleep 2s # Slows down script to make terminal output more readable
# If bioinfo-notebook/ is not in the home directory...
if [ ! -d ~/bioinfo-notebook/ ];
then
	echo ERROR: bioinfo-notebook/ is not in the home directory
	echo The home directory is $HOME
	echo Please move the bioinfo-notebook/ directory to the home directory,
	echo or create a copy of bioinfo-notebook/ in $HOME
	exit 1
fi

echo Updating Linux software...
sleep 2s # Slows down script to make terminal output more readable
sudo apt-get update

echo Downloading Miniconda installation script...
sleep 2s # Slows down script to make terminal output more readable
# If the Linux system is 64-bit...
if [ "$(uname -m)" == "x86_64" ];
then
	# Download the script to install the 64-bit version of miniconda
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# If the Linux system is not 64-bit...
else
	# Download the script to install the 32-bit version of miniconda
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh -O miniconda.sh
fi

echo Installing Miniconda3...
sleep 2s # Slows down script to make terminal output more readable
bash miniconda.sh -b -p $HOME/miniconda

echo Setting up Miniconda3...
sleep 2s # Slows down script to make terminal output more readable
source "$HOME/miniconda/etc/profile.d/conda.sh"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda

echo Displaying information about current conda install...
sleep 2s # Slows down script to make terminal output more readable
conda info -a

echo Creating the bioinfo-notebook virtual environment using conda...
sleep 2s # Slows down script to make terminal output more readable
# If the Linux system is 64-bit...
if [ "$(uname -m)" == "x86_64" ];
then
	# Create the virtual environment using the explicit spec list
	conda create --name bioinfo-notebook --file ~/bioinfo-notebook/envs/bioinfo-notebook.txt
# If the Linux system is not 64-bit...
else
	# Create the virtual environment using an "environment".yml file
	conda env create -f ~/bioinfo-notebook/envs/bioinfo-notebook.txt
fi

echo Script finished
