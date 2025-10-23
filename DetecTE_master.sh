#!/bin/bash

##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                       DetecTE master : a script to manage several                       #############
############                           detecTE annalysis and merge results                           #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################
echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : starting..."
module load python blast/2.14.0 mafft/7.407

##################################################################################################################
##################################################################################################################
############                                       ARGUMENTS                                         #############
##################################################################################################################
##################################################################################################################
# fasta_files should be filled with path to transcriptome files to extract TE from and merge the extracted TE
# assembler_list should contain the assembler corresponging to each fasta file (same order) or one assembler name use to assemble every fasta
# database_list, same with path to protein database file(s)
# prefix_list, same with prefix for each analysis (corresponding to each fasta). If none given, the prefix would be the name of the fasta
# outdir should be the path to the directory with all outputs 

#### Initialisation ####
fasta_files=()
assembler_list=()
database_list=()
prefix_list=()
outdir="./"
dataset_name="Det"

#### usage function ####
usage() {
    echo "Usage: $0 -f <transcrits.fa> -d <db.fa> [OPTIONS]"
    echo "  -f, --fasta <transcrits.fa>        : List of paths to transcriptome files to extract TE from and merge the extracted TE (mandatory)"
    echo "  -d, --database <db.fa>             : Path to the protein database file in fasta format, or list of path is one database for each transcriptome (mandatory)"
    echo "  -o, --outdir <output_dir>          : Path to the output directory (mandatory)"
    echo "  -p, --prefix <prefix>              : Prefix for output files (optionnal)"
    echo "  -n, --dataset_name <dataset_name>  : Prefix merge outputs (optionnal)"
    echo "  -a, --assembler <ASSEMBLER>        : Assembler name in uppar case (TRINITY or SPADES), or list of names if different for the transcriptomes (mandatory)"
    exit 1
}

current_arg=""
#### parsing arguments ####
for arg in "$@"; do
  case "$arg" in
    -f | --fasta)
      current_arg="fasta"
      ;;
    -o | --outdir)
      current_arg="outdir"
      ;;
    -a | --assembler)
      current_arg="assembler"
      ;;
    -d | --database)
      current_arg="database"
      ;;
    -p | --prefix)
      current_arg="prefix"
      ;;
    -n | --dataset_name)
      current_arg="dataset_name"
      ;;
    *)
      if [ "$current_arg" == "fasta" ]; then
          fasta_files+=("$arg")
      elif [ "$current_arg" == "assembler" ]; then
          assembler_list+=("$arg")
      elif [ "$current_arg" == "outdir" ]; then
          outdir=("$arg")
      elif [ "$current_arg" == "database" ]; then
          database_list+=("$arg")
      elif [ "$current_arg" == "prefix" ]; then
          prefix_list+=("$arg")
      elif [ "$current_arg" == "dataset_name" ]; then
          dataset_name="$arg"
      else
          echo "Erreur: Argument inconnu $arg"
          usage
      fi
      ;;
  esac
done

#### Check argument format ####
## fasta 
if [ ${#fasta_files[@]} -eq 0 ]; then # check argument length (here equals 0)
  echo "Error: --fasta must be given with at least one file"
  usage
else
  dataset_size=${#fasta_files[@]}
fi
## assembler
if [ ${#assembler_list[@]} -eq 0 ]; then
  echo "Error: --assembler must be given with at least one file"
  usage
elif [ ${#assembler_list[@]} -eq 1 ]; then # if only one assembler, it is applied to all transcriptomes
  for i in `seq 2 1 ${dataset_size}` ; do
    assembler_list+=("${assembler_list[0]}")
  done
elif [ ! ${#assembler_list[@]} -eq ${dataset_size} ]; then
  echo "Error: assembler list size should be either one or the size of the fasta list"
  echo "Assembler list size : ${#assembler_list[@]}"
  echo "fasta list size : ${#fasta_files[@]}"
  usage
fi
## database
if [ ${#database_list[@]} -eq 0 ]; then
  echo "Erreur: --database must be given with at least one file"
  usage
elif [ ${#database_list[@]} -eq 1 ]; then # if only one database, it is applied to all transcriptomes
  for i in `seq 2 1 ${dataset_size}` ; do
    database_list+=("${database_list[0]}")
  done
elif [ ! ${#database_list[@]} -eq ${dataset_size} ]; then
  echo "Error: database list size should be either one or the size of the fasta list"
  echo "Database list size : ${#database_list[@]}"
  usage
fi
## prefix
if [ ${#prefix_list[@]} -eq 0 ]; then # no prefix given, the name of the transcriptomes files are used as prefix
  for i in `seq 0 1 $((${dataset_size}-1 ))` ; do
    prefix_list+=("`echo "$(basename ${fasta_files[$i]})" | cut -d '.' -f 1`")
  done
elif [ ! ${#prefix_list[@]} -eq ${dataset_size} ]; then
  echo "Error: if given, the size of the prefix list should be the size of the fasta list"
  echo "Prefix list size : ${#prefix_list[@]}"
  usage
fi

echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Parameters have been checked"

#### Show arguments values ####
echo "Queries:"
for file in "${fasta_files[@]}"; do
  echo "  $file"
done
echo "Assemblers:"
for assembler in "${assembler_list[@]}"; do
  echo "  $assembler"
done
echo "Databases:"
for database in "${database_list[@]}"; do
  echo "  $database"
done

echo "Dataset name : ${dataset_name}"
outdir="${outdir%/}_$(date +"%Y-%m-%d")/" # if needed, add "/" behind outdir

##################################################################################################################
##################################################################################################################
############                                        DETECTE                                          #############
##################################################################################################################
##################################################################################################################
echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Submitting detecTE job for each set..."

#### executing detecTE for each transcriptome #### 
# The job ID is stored to wait for the end of the job to continue
job_list=""
for set in `seq 0 1 $((${dataset_size}-1 ))`; do
  mkdir -p ${outdir}${prefix_list[$set]}
  sed "/^#SBATCH -o/s|/shared.*$|${outdir}${prefix_list[$set]}/Report_detecTE.%j.out|" detecTE.sh  > ${outdir}/${prefix_list[$set]}/tmp_detecTE.sh 
  job_detecte=$(sbatch -J ${prefix_list[$set]}_detecTE ${outdir}/${prefix_list[$set]}/tmp_detecTE.sh --fasta "${fasta_files[$set]}" --database "${database_list[$set]}" --assembler "${assembler_list[$set]}" --prefix "${prefix_list[$set]}" --outdir "${outdir}${prefix_list[$set]}"  | awk '{print$4}')
  job_list+="|${job_detecte}"
done
job_list=$(echo ${job_list} | cut -c 2-)
echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Submission done :D"

#### wait for all detecTE jobs to end ####
while [[ $(squeue -u llamothe | grep -E ${job_list} | wc -l) -gt 0 ]]; do
    sleep 10 
done
echo -e "$(date +"%Y-%m-%d %H:%M:%S") : All detecTE done :D"

##################################################################################################################
##################################################################################################################
############                                         MERGE                                           #############
##################################################################################################################
##################################################################################################################
echo -e "\n$(date +"%Y-%m-%d %H:%M:%S") : Merging all results..."

#### grouping the TE extracted from each analysis in one file ####
echo -n > ${outdir}${dataset_name}_TE_before_merging.fasta
for set in `seq 0 1 $((${dataset_size}-1 ))`; do
  cat ${outdir}/${prefix_list[$set]}/${prefix_list[$set]}_confirmedTE.fasta >> ${outdir}${dataset_name}_TE_before_merging.fasta
done
for i in `ls ${outdir}`; do 
  chmod a+rwx ${outdir}$i
done

#### TRIMING TE UNDER 500 BP ####
python TE-trimming.py --fasta ${outdir}${dataset_name}_TE_before_merging.fasta --outfile ${outdir}${dataset_name}_trimmed_TE_before_merging.fasta

echo -e "$(date +"%Y-%m-%d %H:%M:%S") : Blastn all vs. all..."


#### Blasting all the TE vs themself to find the identical TE #### 
mkdir -p ${outdir}tmp_merge
makeblastdb -in ${outdir}${dataset_name}_trimmed_TE_before_merging.fasta -dbtype nucl -out ${outdir}tmp_merge/all_datablastn
blastn -query ${outdir}${dataset_name}_trimmed_TE_before_merging.fasta -db ${outdir}tmp_merge/all_datablastn -outfmt 6 -evalue 1E-21 -out ${outdir}tmp_merge/${dataset_name}_blast_all.out
echo -e "$(date +"%Y-%m-%d %H:%M:%S") : done :D"


#### Extracting connected components from blastn (~ TE famillies) ####
#### writing count_table with number of familly by superfamilly ####
mkdir -p ${outdir}tmp_merge/famillies/alignment
mkdir -p ${outdir}superfamillies
python connected_component.py --fasta ${outdir}${dataset_name}_trimmed_TE_before_merging.fasta --blast ${outdir}tmp_merge/${dataset_name}_blast_all.out  --outdir ${outdir}tmp_merge/famillies/${dataset_name}_ --prefix ${dataset_name} 
mv ${outdir}tmp_merge/famillies/${dataset_name}_count_table.csv ${outdir}
for i in `ls ${outdir}tmp_merge/famillies/`; do 
  chmod a+rwx ${outdir}tmp_merge/famillies/$i
done
echo -e "$(date +"%Y-%m-%d %H:%M:%S") : All detecTE done :D"


#### moving unknown in trash folder ####
mkdir -p ${outdir}tmp_merge/famillies/trash
for i in `ls ${outdir}tmp_merge/famillies/ | grep 'Inter'`; do
  mv ${outdir}tmp_merge/famillies/$i ${outdir}tmp_merge/famillies/trash
done