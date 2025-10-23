#!/bin/bash

##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                          detecTE : a script to detect and count                         #############
############                 the transposons families in an assembled transcriptome                  #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################
module load python blast/2.14.0

##################################################################################################################
############                                       ARGUMENTS                                         #############
##################################################################################################################

#### Default values for arguments ####
fasta=""
outdir="."
prefix=""
database=""
assembler="TRINITY"

#### usage function ####
usage() {
    echo "Usage: $0 -f <transcrits.fa> -d <db.fa> [OPTIONS]"
    echo "  -f, --fasta <transcrits.fa>  : Path to the transcriptome file in fasta format (mandatory)"
    echo "  -d, --database <db.fa>       : Path to the protein database file in fasta format (mandatory)"
    echo "  -o, --outdir <output_dir>    : Path to the output directory (optionnal)"
    echo "  -p, --prefix <prefix>        : Prefix for output files (optionnal)"
    echo "  -a, --assembler <ASSEMBLER>  : Assembler name in uppar case (TRINITY or SPADES)) (optional)"
    exit 1
}

#### Parse command-line options ####
while [[ $# -gt 0 ]]; do
    case $1 in
    -f | --fasta)
        fasta=$2
        shift 2
        ;;
    -d | --database)
        database=$2
        shift 2
        ;;
    -o | --outdir)
        outdir=$2
        shift 2
        ;;
    -p | --prefix)
        prefix=$2
        shift 2
        ;;
    -a | --assembler)
        assembler=$2
        shift 2
        ;;
    -h | --help)
        usage
        ;;
    *)
        echo "Unknown option: $1"
        usage
        ;;
    esac
done

#### Check if required arguments are provided ####
if [ -z "$fasta" ] || [ -z "$database" ]; then
    echo "Error: Missing mandatory argument(s)"
    usage
fi

#### optionnal paramters default value ####
if [ -z "$prefix" ]; then
    prefix=$(basename "$fasta" | awk -F "." '{print $1}')
fi
outdir="${outdir%/}/"

#### Showing parameters value ####
echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Parameters have been checked"
echo "Fasta: $fasta"
echo "Database: $database"
echo "Output Directory: $outdir"
echo "Prefix: $prefix"
echo "Assembler: $assembler"
mkdir -p ${outdir}tmp/

##################################################################################################################
############                                   FILTER FORMAT SPLIT                                   #############
##################################################################################################################

echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : formatting and filtering transcriptome, splitting databse..."

python filter_and_format.py --fasta ${fasta} --outdir ${outdir}tmp/ --assembler ${assembler} --prefix ${prefix}
fasta=${outdir}tmp/${prefix}_all_transcripts.fasta
mkdir -p ${outdir}tmp/TEdb-split
dbname="$(basename ${database} | awk -F '.' '{print $1}')"
seq_number=$(($(grep ">" ${database} | wc -l) / 20 )) 
python split_seqs.py --fasta ${database} --outdir ${outdir}tmp/TEdb-split/ --prefix ${dbname} --seq_num ${seq_number}

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Done :D"

##################################################################################################################
############                                         TBLASTN                                         #############
##################################################################################################################

echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Preparing Tblastn..."

mkdir -p ${outdir}tmp/transcript_index
makeblastdb -in ${fasta} -dbtype nucl -out ${outdir}tmp/transcript_index/${prefix}_datatblastn
querynumber=$(ls ${outdir}tmp/TEdb-split | wc -l) ; querynumber=$((${querynumber} - 1 ))
sed "/^#SBATCH --array=0-0$/s/0-0/0-${querynumber}/" tblastn_det.sh > ${outdir}tmp/tblastn_det_tmp.sh 
mkdir -p ${outdir}tmp/TBLASTNresults/reports
job_tblastn=$(sbatch -J ${prefix}_tblastn ${outdir}tmp/tblastn_det_tmp.sh -d ${outdir}tmp/transcript_index/${prefix}_datatblastn -q ${outdir}tmp/TEdb-split/ -o ${outdir}tmp/TBLASTNresults/${prefix} | cut -d " " -f 4)

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Tblastn running..."

while [[ $(squeue -j $job_tblastn -h | wc -l) -gt 0 ]]; do
    sleep 5 
done

echo -n > ${outdir}tmp/TBLASTNresults/${prefix}_all_TBLASTN.out
for i in `ls ${outdir}tmp/TBLASTNresults/ | grep -v 'reports' | grep -v 'all'` ; do  
    cat ${outdir}tmp/TBLASTNresults/$i >> ${outdir}tmp/TBLASTNresults/${prefix}_all_TBLASTN.out
done

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Done :D"

##################################################################################################################
############                                       FILTER SPLIT                                      #############
##################################################################################################################

echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : second filter tblastn hits, splitting putative TE-transcripts file..."

python second_filter.py --fasta ${fasta} --input ${outdir}tmp/TBLASTNresults/${prefix}_all_TBLASTN.out --prefix ${prefix} --outdir ${outdir}tmp/
mkdir -p ${outdir}tmp/putative-TE-split
seq_number=$(($(grep ">" ${outdir}tmp/${prefix}_putative_TE-transcripts.fasta | wc -l) / 20 )) 
python split_seqs.py --fasta ${outdir}tmp/${prefix}_putative_TE-transcripts.fasta --outdir ${outdir}tmp/putative-TE-split/ --prefix ${prefix}_putative_TE-transcripts --seq_num ${seq_number}

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Done :D"


##################################################################################################################
############                                          BLASTX                                         #############
##################################################################################################################


echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Preparing Blastx..."

mkdir -p ${outdir}tmp/TE_index
makeblastdb -in ${database} -dbtype prot -out ${outdir}tmp/TE_index/$(basename $database | cut -d "." -f 1)_datablastx
querynumber=$(ls ${outdir}tmp/putative-TE-split | wc -l) ; querynumber=$((${querynumber} - 1 ))
sed "/^#SBATCH --array=0-0$/s/0-0/0-${querynumber}/" blastx_det.sh  > ${outdir}tmp/blastx_det_tmp.sh 
mkdir -p ${outdir}tmp/BLASTXresults/reports
job_blastx=$(sbatch -J ${prefix}_blastx ${outdir}tmp/blastx_det_tmp.sh -d ${outdir}tmp/TE_index/$(basename $database | cut -d "." -f 1)_datablastx -q ${outdir}tmp/putative-TE-split/ -o ${outdir}tmp/BLASTXresults/${prefix} | cut -d " " -f 4)

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Blastx running..."

while [[ $(squeue -j $job_blastx -h | wc -l) -gt 0 ]]; do
    sleep 5 
done

echo -n > ${outdir}tmp/${prefix}_all_BLASTX.out
for i in `ls ${outdir}tmp/BLASTXresults/ | grep -v reports` ; do  
    cat ${outdir}tmp/BLASTXresults/$i >> ${outdir}tmp/${prefix}_all_BLASTX.out
done

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Done :D"



##################################################################################################################
############                                 CROSSMATCH & EXTRACTION                                 #############
##################################################################################################################

echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Crossmatch and TE extraction..."

python cross_extract_det.py --fasta ${outdir}tmp/${prefix}_putative_TE-transcripts.fasta  --blastx ${outdir}tmp/${prefix}_all_BLASTX.out --tblastn ${outdir}tmp/${prefix}_TBLASTN_filtered.out --prefix ${prefix} --outdir ${outdir} 

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Done :D"
