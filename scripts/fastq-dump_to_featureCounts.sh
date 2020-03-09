#! /bin/bash

# Help/usage text
usage="$(basename "$0") [options] -a|--annotation <annotation_file> \
	-f|--fasta <fasta_file> <SRR ID(s)> \n
\n
This script downloads FASTQ reads from NCBI's SRA, aligns them to an annotated \n
genome using bowtie2, and generates gene count table(s) using featureCounts.\n
It can take a single SRR ID as an input, or multiple SRR IDs separated by\n
spaces.\n
\n
Required arguments: \n
\t      -a | --annotation\t     input genome annotation file \n
\t      -f | --fasta\t\t        input FASTA file for annotated genome \n
\t      SRR ID(s)\t\t           Sequence Read Archive Run ID(s) (SRR...) \n
\n
Options: \n
\t      -h | --help\t\t         show this help text and exit \n
\t      -p | --processors\t	number (n) of processors to use (default: 1) \n
\t      --fastq-dump\t\t        use 'fastq-dump' instead of the 'fasterq-dump'\n
\t      --verbose\t\t           make output of script more verbose
"

# Setting FASTQDUMP to 0
# This will be changed to "1" if --fastq-dump is given as an argument,
# resulting in fastq-dump being used instead of the default fasterq-dump
FASTQDUMP=0

# Setting VERBOSE to 0
# This will be changed to "1" if --verbose is given as an argument,
# resulting in more verbose script output
VERBOSE=0

# Setting default number of PROCESSORS to use
PROCESSORS=1

# Creating an empty variable for SRRs to be downloaded and aligned to genome
SRRs=""

# Print usage instructions if script is called without any arguments
if [ "$1" = "" ] ; then
  echo -e "ERROR: please provide input files. \n"
  echo -e $usage
  exit 1
fi

# Iterating through the input arguments with a while loop
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit
			;;
		-a|--annotation)
			ANNOTATION=$2
			shift 2
			;;
		-f|--fasta)
			FASTA=$2
			shift 2
			;;
		-p|--processors)
			PROCESSORS=$2
			shift 2
			;;
		--fastq-dump)
			FASTQDUMP=1
			shift
			;;
        --verbose)
            VERBOSE=1
            shift
            ;;
		--) # end argument parsing
			shift
			break
			;;
		-*|--*) # unsupported flags
			echo -e "ERROR: $1 is an invalid option. \n" >&2
			echo -e $usage
			exit 1
			;;
		*) # preserve SRR ID(s) as positional arguments
			SRRs="$SRRs $1"
			shift
			;;
	esac
done

# Beginning the main body of the script
# The sleep commands ("sleep 1s", "sleep 2s") slow down the script to make
# the output more readable in real-time

echo -e		~~~~~~~~~~~~~  F A S T Q - D U M P t o F E A T U R E C O U N T S  ~~~~~~~~~~~~~
echo Script started: $(date)

# Loop through the input SRR IDs
for SRR in $SRRs
do
	printf "\n"
	echo ================================================================================
	echo SRR ID: $SRR
	sleep 1s
	echo Reference genome annotation: $ANNOTATION
	sleep 1s
	echo Reference genome multi-FASTA file: $FASTA
	echo ================================================================================
	sleep 1s
	
	if [ $VERBOSE -eq "1" ]
	then
        	printf "\n"
        	echo Listing files in directory ...
        	sleep 1s
        	ls
        	sleep 2s
	fi


	if [ $FASTQDUMP -eq "1" ]
	then
        	if [ $VERBOSE -eq "1" ]
		then
            		echo Downloading compressed FASTQ reads using fastq-dump...
        	fi
		until fastq-dump --gzip --skip-technical --readids --read-filter pass \
		--dumpbase --split-3 --clip $SRR; do
			echo fastq-dump failed, retrying in 10 seconds...
	    		sleep 10s
		done
	else
        	if [ $VERBOSE -eq "1" ]
		then
            		echo Downloading FASTQ reads using fasterq-dump...
        	fi
		until fasterq-dump --progress --threads $PROCESSORS $SRR; do
			echo fasterq-dump failed, retrying in 10 seconds...
			rm -r fasterq.tmp.*
			sleep 10s
		done
	fi

	if [ $VERBOSE -eq "1" ]
	then
        	sleep 1s
        	echo Listing files in directory after downloading reads...
        	sleep 1s
        	ls
        	sleep 2s
    	fi

	# Checking if bowtie2 index of FASTA file exists before creating bowtie2 index
	# If bowtie2_$FASTA.1.bt2 (one of the bowtie2 index files) does not exist...
	if [ ! -f bowtie2_$FASTA.1.bt2 ]
	# ...then create the bowtie2_$FASTA index
	then
        	if [ $VERBOSE -eq "1" ]
		then
            		echo Indexing reference genome FASTA file using bowtie2-build...
			sleep 2s
		fi
	    	bowtie2-build $FASTA bowtie2_$FASTA
	    	if [ $VERBOSE -eq "1" ]
		then
            		sleep 1s
            		echo Listing files in directory after running bowtie2-build...
            		sleep 1s
            		ls
            		sleep 2s
        	fi
	# Otherwise, print a message confirming that it exists
	else
        	if [ $VERBOSE -eq "1" ]
		then
            		echo The bowtie2 index bowtie2_$FASTA exists
            		sleep 1s
	    	fi
	fi

	if [ $VERBOSE -eq "1" ]
	then
        	echo Aligning reads to reference genome using bowtie2...
        	sleep 2s
    	fi

	# Checking if fastq-dump or fasterq-dump was used, as this will result
	# in different filenames
	if [ $FASTQDUMP -eq "1" ]
	then
		bowtie2 -p $PROCESSORS --no-unal -x bowtie2_$FASTA \
		-1 $SRR\_pass_1.fastq.gz -2 $SRR\_pass_2.fastq.gz \
		-S $SRR\_$FASTA.sam
	else
		bowtie2 -p $PROCESSORS --no-unal -x bowtie2_$FASTA \
		-1 $SRR\_1.fastq -2 $SRR\_2.fastq \
		-S $SRR\_$FASTA.sam
	fi

	if [ $VERBOSE -eq "1" ]
	then
		sleep 1s
        	echo Listing files in directory after running bowtie2...
        	sleep 1s
        	ls
        	sleep 2s

        	echo Converting alignment from SAM to BAM format using samtools view...
        	sleep 2s
    	fi
	samtools view -@ $PROCESSORS -Sb $SRR\_$FASTA.sam \
	> $SRR\_$FASTA.bam

	if [ $VERBOSE -eq "1" ]
	then
        	sleep 1s
	        echo Listing files in directory after running samtools view...
        	sleep 1s
        	ls
        	sleep 2s

        	echo Sorting the BAM file using samtools sort...
        	sleep 2s
    	fi
	samtools sort -@ $PROCESSORS $SRR\_$FASTA.bam \
	-o sorted_$SRR\_$FASTA.bam

	if [ $VERBOSE -eq "1" ]
	then
        	sleep 1s
        	echo Listing files in directory after running samtools sort...
        	sleep 1s
        	ls
        	sleep 2s
    
        	echo Generating count table using featureCounts...
        	sleep 2s
    	fi
	featureCounts -p -s 2 -T $PROCESSORS -a $ANNOTATION \
	-o feature_counts_$SRR\_$FASTA.txt \
	sorted_$SRR\_$FASTA.bam

	if [ $VERBOSE -eq "1" ]
	then
        	sleep 1s
        	echo Listing files in directory after running featureCounts...
        	sleep 1s
        	ls
        	sleep 2s
	
        	echo Results written to feature_counts_$SRR\_$FASTA.txt
        	sleep 2s

        	echo Head of feature_counts_$SRR\_$FASTA.txt
        	sleep 2s
        	head feature_counts_$SRR\_$FASTA.txt
        	sleep 2s

        	echo Tail of feature_counts_$SRR\_$FASTA.txt
        	sleep 2s
        	tail feature_counts_$SRR\_$FASTA.txt
        	sleep 2s
    	fi
done

echo Script finished: $(date)
