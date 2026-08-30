#!/bin/bash

### Default values for arguments ###
querydir=""
outdir="."
prefix="res"
database=""

### usage function ###
usage() {
    echo "Usage: $0 -f <transcrits.fa> -d <db.fa> [OPTIONS]"
    echo "  -q, --querydir <file.fa>     : Path to the directory with a list of protein query files in fasta format (mandatory)"
    echo "  -d, --database <db.fa>       : Path to the indexed nucleotides database file in fasta format (mandatory)"
    echo "  -p, --prefix <string>        : Prefix of the output files (optionnal)"
    echo "  -o, --outfile <output_dir>   : Path to the output file (mandatory)"
    exit 0
}

### Parse command-line options ###
while [[ $# -gt 0 ]]; do
    case $1 in
    -q | --query)
        querydir=$2
        shift 2
        ;;
    -d | --database)
        database=$2
        shift 2
        ;;
    -o | --outfile)
        outdir=$2
        shift 2
        ;;
    -p | --prefix)
    outdir=$2
    shift 2
    ;;
    -h | --help)
        usage
        shift 2
        ;;
    *) # unknown option
        echo "Unknown option: $key"
        usage
        shift 2
        ;;
    esac
done

### Check if required arguments are provided ###
if [ -z "$querydir" ] || [ -z "$database" ] || [ -z "$outdir" ]; then
    echo "Error: Missing mandatory argument(s)"
    usage
fi


### Set the files needed by blast ###
fastaname=${querydir}
outfilename="${outdir}${prefix}_tblastn-$SLURM_ARRAY_TASK_ID.out"
echo "-db : ${database} "
echo "-query : ${fastaname[$SLURM_ARRAY_TASK_ID]}  "
echo "-out : ${outfilename} "
### Run blast ###
tblastn -db ${database} -query ${fastaname[$SLURM_ARRAY_TASK_ID]} -outfmt 6 -evalue 1E-21 -out ${outfilename}
# tblastn -db ${database} -query ${fastaname[$SLURM_ARRAY_TASK_ID]} -outfmt 6 -evalue 1E-10 -out ${outfilename}


echo -e "\n\n$(date +"%Y-%m-%d %H:%M:%S") : i'm done"