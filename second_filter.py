##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                                      Sequence size filter :                             #############
############                   a script to filter the too sequences in a fasta file                  #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import numpy as np
import argparse
parser = argparse.ArgumentParser(prog="transcript_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--input", "-i", required=True, type=str, action="store", help="blast output (fmt 6)")
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help="pPrefix before the output files")
args = parser.parse_args()


# Removes the false positive. The sequence is kept if :
#       - percentage of identity > 25%
#       - alignement length > 50 aa
#       - number of gaps < 15
# more about output format : https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# output :
#       - filtered_seqs = list of ID of sequences meeting all criteria
#       - [prefix]_BLASTX_filtered.out = same as [prefix]_BLASTX.out without the sequence deleted by the filter
#       - [prefix]_TEputatifs.fasta = .fasta file with filtered sequences

infile = open(args.input, "r")
outfile = open("{}{}_TBLASTN_filtered.out".format(args.outdir, args.prefix), "w")
filtered_seqs = []
for line in infile:
    item = line[:-1].split("\t")
    title = "_".join(line[:-1].split("\t")[1].split("_"))
    if (float(item[2]) > 25 and float(item[3]) > 50 and float(item[5]) < 15):
        filtered_seqs.append(title)
        outfile.write(line)
infile.close()
outfile.close()
filtered_seqs = np.unique(filtered_seqs)


#### CREATING FASTA WITH PUTATIVE TEs ####
infile = open(args.fasta, "r")
outfile = open("{}{}_putative_TE-transcripts.fasta".format(args.outdir, args.prefix), "w")
found = False
for line in infile:
    if line[0] == ">":
        line="_".join(line[1:].split("_")[-8:])
        if line[:-1] in filtered_seqs:
            found = True
            titre = line.split("_")
            line = ">" + titre[-2] + "_" + titre[-1]
        else:
            found = False
    if found:
        outfile.write(line)

infile.close()
outfile.close()



