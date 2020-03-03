#! /bin/bash

# Defining script arguments using getopts...
while getopts :i:a:f:p: OPT; do
    case "$OPT" in
        i) SRR=$OPTARG ;;
        a) ANNOTATION=$OPTARG ;;
        f) FASTA=$OPTARG ;;
        p) PROCESSORS=$OPTARG ;;
    esac
done

# Print usage instructions if script is called without any arguments
if [ "$1" = "" ] ; then
  echo "-i  Sequence Read Archive Run ID (SRR...)"
  echo "-a  Reference genome annotation file"
  echo "-f  Reference genome multi-FASTA file"
  echo "-p  Number of processors/threads to use"
  exit 0
fi

printf "\n"
echo ===========================================================================
echo fastq-dump_to_featureCounts.sh
echo Script started: $(date)
sleep 1s
echo SRR ID: $SRR
sleep 1s
echo Reference genome annotation: $ANNOTATION
sleep 1s
echo Reference genome multi-FASTA file: $FASTA
echo ===========================================================================
sleep 1s

printf "\n"
echo Listing files in directory ...
ls
sleep 2s

echo Downloading compressed FASTQ reads using fastq-dump... ~~~~~~~~~~~~~~~~~~~~
until fastq-dump --gzip --skip-technical --readids --read-filter pass \
--dumpbase --split-3 --clip $SRR; do
    echo fastq-dump failed, retrying in 10 seconds...
    sleep 10s
done

sleep 2s

echo Indexing reference genome FASTA file using bowtie2-build ~~~~~~~~~~~~~~~~~~
sleep 2s
bowtie2-build $FASTA bowtie2_$FASTA

sleep 2s

echo Aligning reads to reference genome using bowtie2 ~~~~~~~~~~~~~~~~~~~~~~~~~~
sleep 2s
bowtie2 -p $PROCESSORS --no-unal -x bowtie2_$FASTA \
-1 $SRR\_pass_1.fastq.gz -2 $SRR\_pass_2.fastq.gz \
-S $SRR\_$FASTA.sam

sleep 2s

echo Converting alignment from SAM to BAM format using samtools view ~~~~~~~~~~~
sleep 2s
samtools view -@ $PROCESSORS -Sb $SRR\_$FASTA.sam \
> $SRR\_$FASTA.bam

sleep 2s

echo Sorting the BAM file using samtools sort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sleep 2s
samtools sort $SRR\_$FASTA.bam \
-o sorted_$SRR\_$FASTA.bam

sleep 2s

echo Generating count table using featureCounts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sleep 2s
featureCounts -p -s 2 -T $PROCESSORS -a $ANNOTATION \
-o feature_counts_$SRR\_$FASTA.txt \
sorted_$SRR\_$FASTA.bam

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

echo Script finished: $(date)
sleep 2s

echo Exiting shell...
sleep 2s
exit
