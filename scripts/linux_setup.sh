#! /bin/bash/

# Help/usage text
usage="$(basename "$0") \n
\n
This script downloads and installs Miniconda3, and uses conda to install \n
the 'bioinfo-notebook' virtual environment. \n
\n
This script requires admin rights, please run as sudo, for example... \n
\t \$ sudo bash $0 \n
\n
Optional arguments: \n
\t      -h | --help\t\t         show this help text and exit \n
"

# Iterating through the input arguments with a while loop
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit 0
			;;
	esac
done

# Changing directory to the home directory ("~" or "$HOME")
cd ~

echo Checking if the bioinfo-notebook environment is already installed...
sleep 2s # Slows down script to make terminal output more readable
if [ -d ~/miniconda/envs/bioinfo-notebook ]; then 
	echo The bioinfo-notebook environment already exists, exiting script.
	exit 0
fi

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

echo Miniconda3 installed, removing installation script...
rm -f miniconda.sh

echo Setting up Miniconda3...
sleep 2s # Slows down script to make terminal output more readable
source "$HOME/miniconda/etc/profile.d/conda.sh"
hash -r
conda config --set always_yes yes --set changeps1 yes
conda update -q conda
conda init

echo Displaying information about current conda installation...
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
	conda env create -f ~/bioinfo-notebook/envs/bioinfo-notebook.yml
fi

echo -e Script finished! \n

echo -e Please restart your Linux system for these changes to take effect. \n

echo The bioinfo-notebook environment can be activated using the command...
echo -e	\$ conda activate bioinfo-notebook
echo A conda virtual environment can be deactivated using the command...
echo -e	\$ conda deactivate
