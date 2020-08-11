#! /bin/bash

# Help/usage text
usage="$(basename "$0") [-h|--help] [-l|--log -p|--processors n] \n
\n
This script downloads FASTQ sequencing reads, aligns them to a reference genome, \n
and creates a VCF file based on this alignment. It is written to be run \n
from the 'bioinfo-notebook/' directory. \n
\n
where: \n
    -h | --help          show this help text and exit \n
    -l | --log           redirect terminal output to a log file \n
    -p | --processors    optional: set the number (n) of processors to use \n
                         (default: 1) \n
"

MAKELOG=false
PROCESSORS=1

# Iterating through the input arguments with a while loop
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo -e $usage
            exit
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
    # CREATING LOG FILE
    # Terminal output directed to the file 'snp_calling_[date]_[time].log'
    exec 3>&1 4>&2
    trap 'exec 2>&4 1>&3' 0 1 2 3
    exec 1>../results/snp_calling_$(date +%Y%m%d_%H%M).log 2>&1
fi

echo "$(date +%Y/%m/%d\ %H:%M) Beginning SNP calling script."

echo Downloading reads...

until fastq-dump --gzip --skip-technical --readids --read-filter pass \
--dumpbase --split-files --clip DRR237290; do
    echo fastq-dump failed, retrying in 10 seconds...
    sleep 10s
done

echo Downloading reference sequence...

curl -s --remote-name --remote-time ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA\
/003/086/655/GCA_003086655.1_ASM308665v1\
/GCA_003086655.1_ASM308665v1_genomic.fna.gz

echo Decompressing reference sequence...

gunzip GCA_003086655.1_ASM308665v1_genomic.fna.gz

echo Indexing reference sequence for bowtie2...

bowtie2-build GCA_003086655.1_ASM308665v1_genomic.fna S_cere_ref_seq

echo Aligning reads to the reference genome...

bowtie2 --no-unal -p $PROCESSORS -x S_cere_ref_seq -1 DRR237290_pass_1.fastq.gz \
	-2 DRR237290_pass_2.fastq.gz -S S_cere_DRR237290_alignment.sam

echo Converting SAM alignment to sorted BAM alignment...

samtools view -@ $PROCESSORS -Sb -o S_cere_DRR237290_alignment.bam \
	S_cere_DRR237290_alignment.sam 

samtools sort -@ $PROCESSORS -O bam -o sorted_S_cere_DRR237290_alignment.bam \
	S_cere_DRR237290_alignment.bam

echo Indexing reference sequence for SAMtools...

samtools faidx GCA_003086655.1_ASM308665v1_genomic.fna

echo Generating genotype variant likelihoods...

samtools mpileup -g -f GCA_003086655.1_ASM308665v1_genomic.fna \
	-o S_cere_DRR237290_full.bcf sorted_S_cere_DRR237290_alignment.bam

echo SNP calling with BCFtools...

bcftools call -O b --threads $PROCESSORS -vc S_cere_DRR237290_full.bcf \
	> S_cere_DRR237290_var.bcf

echo SNP filtering with vcfutils.pl...

bcftools view S_cere_DRR237290_var.bcf | vcfutils.pl varFilter - \
	> S_cere_DRR237290_SNP.vcf

echo Head of VCF file...

head S_cere_DRR237290_SNP.vcf

echo Tail of VCF file...

tail S_cere_DRR237290_SNP.vcf

echo "$(date +%Y/%m/%d\ %H:%M) Script finished."
