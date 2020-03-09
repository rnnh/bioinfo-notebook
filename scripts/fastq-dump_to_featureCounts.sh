#! /bin/bash

# Help/usage text
usage="$(basename "$0") [-h|--help] [-a|--annotation *.gtf -f|--fasta *.fasta -p|--processors n] <SRR ID(s)> \n
\n
This script downloads FASTQ reads from NCBI's SRA, aligns them to an annotated genome using bowtie2, and generates gene count table(s) using featureCounts. It can take a single SRR ID as an input, or multiple SRR IDs separated by spaces. \n
\n
where: \n
    -h | --help          show this help text and exit \n
    -a | --annotation    input genome annotation file \n
    -f | --fasta        input FASTA file for annotated genome \n
    -p | --processors    number (n) of processors to use (optional, default: 1) \n
    SRR ID(s)            Sequence Read Archive Run ID(s) (SRR...) \n
"

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

echo        ~~~ F A S T Q - D U M P t o F E A T U R E C O U N T S ~~~
echo Script started: $(date)

# Loop through the input SRR IDs
for SRR in $SRRs
do
	printf "\n"
	echo ===========================================================================
	echo SRR ID: $SRR
	sleep 1s
	echo Reference genome annotation: $ANNOTATION
	sleep 1s
	echo Reference genome multi-FASTA file: $FASTA
	echo ===========================================================================
	sleep 1s

	printf "\n"
	echo Listing files in directory ...
	sleep 1s
	ls
	sleep 2s

	echo Downloading compressed FASTQ reads using fastq-dump... ~~~~~~~~~~~~~~~~~~~~
	until fastq-dump --gzip --skip-technical --readids --read-filter pass \
	--dumpbase --split-3 --clip $SRR; do
	    echo fastq-dump failed, retrying in 10 seconds...
	    sleep 10s
	done

	sleep 1s
	echo Listing files in directory after running fastq-dump...
	sleep 1s
	ls
	sleep 2s

	# Checking if bowtie2 index of FASTA file exists before creating bowtie2 index
	# If bowtie2_$FASTA.1.bt2 (one of the bowtie2 index files) does not exist...
	if [ ! -f bowtie2_$FASTA.1.bt2 ]
	# ...then create the bowtie2_$FASTA index
	then
	    echo Indexing reference genome FASTA file using bowtie2-build ~~~~~~~~~~~~~~~~~~
	    sleep 2s
	    bowtie2-build $FASTA bowtie2_$FASTA
	    sleep 1s
	    echo Listing files in directory after running bowtie2-build...
	    sleep 1s
	    ls
	    sleep 2s
	# Otherwise, print a message confirming that it exists
	else
	    echo The bowtie2 index bowtie2_$FASTA exists
	    sleep 1s
	fi

	echo Aligning reads to reference genome using bowtie2 ~~~~~~~~~~~~~~~~~~~~~~~~~~
	sleep 2s
	bowtie2 -p $PROCESSORS --no-unal -x bowtie2_$FASTA \
	-1 $SRR\_pass_1.fastq.gz -2 $SRR\_pass_2.fastq.gz \
	-S $SRR\_$FASTA.sam

	sleep 1s
	echo Listing files in directory after running bowtie2...
	sleep 1s
	ls
	sleep 2s

	echo Converting alignment from SAM to BAM format using samtools view ~~~~~~~~~~~
	sleep 2s
	samtools view -@ $PROCESSORS -Sb $SRR\_$FASTA.sam \
	> $SRR\_$FASTA.bam

	sleep 1s
	echo Listing files in directory after running samtools view...
	sleep 1s
	ls
	sleep 2s

	echo Sorting the BAM file using samtools sort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	sleep 2s
	samtools sort -@ $PROCESSORS $SRR\_$FASTA.bam \
	-o sorted_$SRR\_$FASTA.bam

	sleep 1s
	echo Listing files in directory after running samtools sort...
	sleep 1s
	ls
	sleep 2s

	echo Generating count table using featureCounts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	sleep 2s
	featureCounts -p -s 2 -T $PROCESSORS -a $ANNOTATION \
	-o feature_counts_$SRR\_$FASTA.txt \
	sorted_$SRR\_$FASTA.bam

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
done

echo Script finished: $(date)
