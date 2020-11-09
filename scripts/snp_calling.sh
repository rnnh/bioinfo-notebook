#! /bin/bash

# Help/usage text
usage="$(basename "$0") [-h|--help] [-1|--one -2|--two -r|--reference] \n
[-d|--demo] [-l|--log -p|--processors n] \n
\n
This script aligns sequencing reads to a reference genome, and finds genetic \n
variants (SNPs/indels) based on this alignment, which are written to a variant\n
call format (VCF) file. This script should be called from the \n
'bioinfo-notebook/' directory.\n
\n
Calling this script with the argument '-d' or '--demo' will run this script \n
using Saccharomyces cerevisiae FASTQ sequencing reads and a Saccharomyces \n
cerevisiae reference genome, which will be downloaded from NCBI. \n
\n
The programs required for this script are in the 'bioinfo-notebook' conda \n
environment (bioinfo-notebook/envs/bioinfo-notebook.yml or \n
bioinfo-notebook/envs/bioinfo-notebook.txt). \n
\n
arguments: \n
\t  -h | --help\t\t          show this help text and exit \n
\t  -1 | --one\t\t           forward reads to align with reference sequence \n
\t\t\t\t                     (FASTQ: .fastq or .fastq.gz) \n
\t  -2 | --two\t\t           reverse reads to align with reference sequence \n
\t\t\t\t                     (FASTQ: .fastq or .fastq.gz) \n
\t  -r | --reference\t       reference sequence to align reads against \n
\t\t\t\t                     (FASTA nucleotide file: .fna) \n
\t  -d | --demo\t\t          run the script with demonstration inputs\n
\t  -l | --log\t\t           redirect terminal output to a log file in the \n
\t\t\t\t                     directory bioinfo-notebook/results/ \n
\t  -p | --processors\t      optional: set the number (n) of processors to \n
\t\t\t\t                     use (default: 1) \n
"

MAKELOG=false
PROCESSORS=1
DEMO=false
ONE=""
TWO=""
REFERENCE=""

# Iterating through the input arguments with a while loop
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo -e $usage
            exit
            ;;
        -1|--one)
            ONE=$2
            shift 2
            ;;
        -2|--two)
            TWO=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -d|--demo)
            DEMO=true
            shift 1
            ;;
        -l|--log)
            MAKELOG=true
            shift 1
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
    esac
done

cd data

if $MAKELOG ; then
    # Creating results directory, if it does not already exist
    if [ ! -d ../results ]; then
        mkdir ../results
    fi
    # CREATING LOG FILE
    # Terminal output directed to the file 'snp_calling_[date]_[time].log'
    exec 3>&1 4>&2
    trap 'exec 2>&4 1>&3' 0 1 2 3
    exec 1>../results/snp_calling_$(date +%Y%m%d_%H%M).log 2>&1
fi

echo "$(date +%Y/%m/%d\ %H:%M) Beginning SNP calling script."

if $DEMO ; then
    echo Downloading reads...
    until fastq-dump --gzip --skip-technical --readids --read-filter pass \
    --dumpbase --split-files --clip DRR237290; do
        echo fastq-dump failed, retrying in 10 seconds...
        sleep 10s
    done

    echo Downloading reference sequence...
    curl -s --remote-name --remote-time ftp://ftp.ncbi.nlm.nih.gov/genomes/all/\
    GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

    echo Decompressing reference sequence...
    gunzip GCF_000146045.2_R64_genomic.fna.gz

    echo Indexing reference sequence for bowtie2...
    bowtie2-build GCF_000146045.2_R64_genomic.fna S_cere_ref_seq

    echo Aligning reads to the reference genome...
    bowtie2 --no-unal -p $PROCESSORS -x S_cere_ref_seq \
    -1 DRR237290_pass_1.fastq.gz -2 DRR237290_pass_2.fastq.gz \
    -S S_cere_DRR237290_alignment.sam

    echo Converting SAM alignment to sorted BAM alignment...
    samtools view -@ $PROCESSORS -Sb \
    -o S_cere_DRR237290_alignment_unsorted.bam S_cere_DRR237290_alignment.sam 

    samtools sort -@ $PROCESSORS -O bam -l 9 -o S_cere_DRR237290_alignment.bam \
    S_cere_DRR237290_alignment_unsorted.bam

    echo Removing redundant alignment files...
    rm -v S_cere_DRR237290_alignment.sam S_cere_DRR237290_alignment_unsorted.bam

    echo Indexing reference sequence for SAMtools...
    samtools faidx GCF_000146045.2_R64_genomic.fna

    echo Generating genotype variant likelihoods with BCFtools...
    bcftools mpileup --max-depth 10000 --threads $PROCESSORS \
    -f GCF_000146045.2_R64_genomic.fna \
    -o S_cere_DRR237290_full.bcf S_cere_DRR237290_alignment.bam

    echo Variant calling with BCFtools...
    bcftools call -O b --threads $PROCESSORS -vc --ploidy 1 -p 0.05 \
    -o S_cere_DRR237290_var_unfiltered.bcf S_cere_DRR237290_full.bcf

    echo Removing redundant BCF file...
    rm -v S_cere_DRR237290_full.bcf

    echo Variant filtering with BCFtools filter...
    bcftools filter --threads $PROCESSORS -i '%QUAL>=20' -O v \
    -o S_cere_DRR237290_var.vcf S_cere_DRR237290_var_unfiltered.bcf

    echo Head of VCF file...
    head S_cere_DRR237290_var.vcf

    echo Tail of VCF file...
    tail S_cere_DRR237290_var.vcf

    echo "$(date +%Y/%m/%d\ %H:%M) Script finished."
fi

echo Indexing reference sequence for bowtie2...
bowtie2-build $REFERENCE reference_seq

echo Aligning reads to the reference genome...
bowtie2 --no-unal -p $PROCESSORS -x reference_seq \
-1 $ONE -2 $TWO -S reference_seq_vs_reads_alignment.sam

echo Converting SAM alignment to sorted BAM alignment...
samtools view -@ $PROCESSORS -Sb \
-o reference_seq_vs_reads_alignment_unsorted.bam \
reference_seq_vs_reads_alignment.sam 

samtools sort -@ $PROCESSORS -O bam -l 9 \
-o reference_seq_vs_reads_alignment.bam \
reference_seq_vs_reads_alignment_unsorted.bam

echo Removing redundant alignment files...
rm -v reference_seq_vs_reads_alignment.sam \
reference_seq_vs_reads_alignment_unsorted.bam

echo Indexing reference sequence for SAMtools...
samtools faidx $REFERENCE

echo Generating genotype variant likelihoods with BCFtools...
bcftools mpileup --max-depth 10000 --threads $PROCESSORS \
-f $REFERENCE -o reference_seq_vs_reads_full.bcf \
reference_seq_vs_reads_alignment.bam

echo Variant calling with BCFtools...
bcftools call -O b --threads $PROCESSORS -vc --ploidy 1 -p 0.05 \
-o reference_seq_vs_reads_var_unfiltered.bcf reference_seq_vs_reads_full.bcf

echo Removing redundant BCF file...
rm reference_seq_vs_reads_full.bcf

echo Variant filtering with BCFtools filter...
bcftools filter --threads $PROCESSORS -i '%QUAL>=20' -O v \
-o reference_seq_vs_reads_var.vcf reference_seq_vs_reads_var_unfiltered.bcf

echo Head of VCF file...
head reference_seq_vs_reads_var.vcf

echo Tail of VCF file...
tail reference_seq_vs_reads_var.vcf

echo "$(date +%Y/%m/%d\ %H:%M) Script finished."
