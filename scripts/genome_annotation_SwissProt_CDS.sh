#! /bin/bash
# https://github.com/rnnh/bioinfo-notebook.git

# Help/usage text
usage="$(basename "$0") [-h|--help] [-d|--demo] [-i|--input] \n
[-l|--log -p|--processors n -e|--email] \n
\n
A script to annotate proteins in a genome assembly, using BLASTx with\n
UniProtKB/Swiss-Prot.\n
\n
When run with the arugment '-d' or '--demo' this script...\n
\n
\t 1. Downloads a Saccharomyces cerevisiae S288C genome assembly, and \n
\t the UniProtKB/Swiss-Prot amino acid sequences. \n
\t 2. Creates a BLAST database from the downloaded Swiss-Prot sequences,\n
\t and searches the S. cerevisiae genome against it using BLASTx with an\n
\t E-value threshold of 1e-100. \n
\t 3. Filters the BLASTx results, removing results with less than 90%\n
\t identity.\n
\t 4. Creates a genome annotation GFF file from these BLASTx results.\n
\t 5. Adds information to the genome annotation from UniProt (protein\n
\t names, KeGG ortholog information, EC numbers, etc.) \n
\n
The end result ('S_cere.gff') is an annotation of the coding sequences (CDS) \n
in the S. cerevisiae genome that are described in UniProtKB/Swiss-Prot. \n
\n
This script can also be run with the argument '-i' or '--input', which is used\n
to specify a FASTA nucleotide file (.fasta or .fna) to annotate, instead of\n
the demo sequence. The end result is also an annotation of the CDS in the input\n
sequence based on UniProtKB/Swiss-Prot, called '<input>.gff'.\n
\n
This script should be called from the 'bioinfo-notebook/' directory.The \n
programs required for this script are in the 'bioinfo-notebook' conda \n
environment (bioinfo-notebook/envs/bioinfo-notebook.yml or \n
bioinfo-notebook/envs/bioinfo-notebook.txt). \n
If the input file is not in the 'bioinfo-notebook/data/' directory, the full \n
file path should be given.\n
\n
arguments: \n
\t  -h | --help\t\t          show this help text and exit \n
\t  -i | --input\t\t         name of input FASTA nucleotide file to annotate \n
\t  -d | --demo\t\t          run the script with demonstration inputs\n
\n
optional arguments:\n
\t  -l | --log\t\t           redirect terminal output to a log file \n
\t  -p | --processors\t      set the number (n) of processors to use\n
\t\t\t\t                     (default: 1) \n
\t  -e | --email\t\t         contact email for UniProt queries
"

MAKELOG=false
PROCESSORS=1
EMAIL="none"
DEMO=false
INPUT=""

# Iterating through the input arguments with a while loop
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo -e $usage
            exit
            ;;
        -i|--input)
            INPUT=$2
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
        -e|--email)
            EMAIL=$2
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
    # Terminal output directed to the file 'genome_annotation_[date]_[time].log'
    exec 3>&1 4>&2
    trap 'exec 2>&4 1>&3' 0 1 2 3
    exec 1>../results/genome_annotation_$(date +%Y%m%d_%H%M).log 2>&1
fi

echo "$(date +%Y/%m/%d\ %H:%M) Beginning genome annotation script."

if $DEMO ; then
    echo Downloading genome FASTA file...
    curl -s -o S_cere.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146\
/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

    echo Decompressing genome FASTA file...
    gunzip S_cere.fna.gz

fi

echo Downloading Swiss-Prot sequences...
curl -s -o uniprot_sprot.fasta.gz ftp://ftp.uniprot.org/pub/databases/uniprot/\
current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | xargs -n 1 \
-P $PROCESSORS

echo Decompressing Swiss-Prot sequences...
gunzip uniprot_sprot.fasta.gz

echo Creating BLAST database...
makeblastdb -dbtype prot -in uniprot_sprot.fasta -out SwissProt

echo Removing Swiss-Prot sequences...
rm -v uniprot_sprot.fasta

if $DEMO ; then
    echo Searching genome FASTA file against Swiss-Prot with BLASTx...
    blastx -num_threads $PROCESSORS -evalue 1e-100 -query S_cere.fna \
    -db SwissProt -outfmt 6 -out blastx_SwissProt_S_cere_unfiltered.tsv

    echo Removing Swiss-Prot database...
    rm -v SwissProt*

echo Filtering BLASTx results with percentage identity less than 90% with awk...
    awk '{ if ($3 >= 90) { print } }' blastx_SwissProt_S_cere_unfiltered.tsv \
    > blastx_SwissProt_S_cere.tsv

    echo Removing unfiltered BLASTx results...
    rm -v blastx_SwissProt_S_cere_unfiltered.tsv

    echo Creating genome annotation GFF file from BLASTx results...
    blast2gff uniprot --fasta-file S_cere.fna blastx_SwissProt_S_cere.tsv \
    S_cere_without_UniProt_info.gff

    echo Adding information to genome annotation from UniProt...
    until add-gff-info uniprot --email $EMAIL --protein-names --enzymes \
    --kegg_orthologs --eggnog --taxon-id S_cere_without_UniProt_info.gff \
    S_cere.gff; do
        echo add-gff-info failed, retrying in 10 seconds...
        rm -v S_cere.gff
        sleep 10s
    done

    echo Removing copy of genome annotation without added UniProt info...
    rm -v S_cere_without_UniProt_info.gff

    echo First line of finished genome annotation...
    head -n 1 S_cere.gff
fi

if [ ! -z $INPUT ]; then
    echo Searching genome FASTA file against Swiss-Prot with BLASTx...
    blastx -num_threads $PROCESSORS -evalue 1e-100 -query $INPUT \
    -db SwissProt -outfmt 6 -out blastx_SwissProt_$INPUT\_unfiltered.tsv

    echo Removing Swiss-Prot database...
    rm -v SwissProt*

echo Filtering BLASTx results with percentage identity less than 90% with awk...
    awk '{ if ($3 >= 90) { print } }' blastx_SwissProt_$INPUT\_unfiltered.tsv \
    > blastx_SwissProt_$INPUT\.tsv

    echo Removing unfiltered BLASTx results...
    rm -v blastx_SwissProt_$INPUT\_unfiltered.tsv

    echo Creating genome annotation GFF file from BLASTx results...
    blast2gff uniprot --fasta-file $INPUT blastx_SwissProt_$INPUT\.tsv \
    $INPUT\_without_UniProt_info.gff

    echo Adding information to genome annotation from UniProt...
    until add-gff-info uniprot --email $EMAIL --protein-names --enzymes \
    --kegg_orthologs --eggnog --taxon-id $INPUT\_without_UniProt_info.gff \
    $INPUT.gff; do
        echo add-gff-info failed, retrying in 10 seconds...
        rm -v $INPUT.gff
        sleep 10s
    done

    echo Removing copy of genome annotation without added UniProt info...
    rm -v $INPUT\_without_UniProt_info.gff

    echo First line of finished genome annotation...
    head -n 1 $INPUT.gff
fi

echo "$(date +%Y/%m/%d\ %H:%M) Script finished."
