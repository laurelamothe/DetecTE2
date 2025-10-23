##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                                      Sequence size filter :                             #############
############                  a script to filter the too long sequences in a fasta file              #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
parser = argparse.ArgumentParser(prog="transcript_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help='Prefix for the output files. If not specified it will take the letters before "." in fasta file name')
parser.add_argument("--assembler", "-a", required=True, type=str, action="store", help='string with the name of the assembler used to generate the transcriptome (TRINITY or SPADES)')
args = parser.parse_args()

infile = open(args.fasta, "r")
outfile = open("{}{}_all_transcripts.fasta".format(args.outdir, args.prefix), "w")
title=""          
seq=""  
length_seq = 0  
i=0      
for line in infile:
    if line[0] == '>':
        i+=1
        if length_seq >= 500:
            outfile.write(title)
            outfile.write(seq)
        if args.assembler == "TRINITY": 
            title = ">" + "-".join(line.split("_")[1:3]) + "_" + "-".join(line.split("_")[3:5]).split(" ")[0] + "\n"
        elif args.assembler == "SPADES":
            title = ">" + "_".join(line.split("_")[-2:])  
        else:
            title =line[1:].split(" ")[0]
            title =  ">" + title.replace("#", "").replace("_", "").replace("/", "-") + "_" + str(i) + "\n"
        seq=""
        length_seq = 0
    else:
        if args.assembler == "TRINITY": 
            seq = seq + "\n".join([line[i:i+60] for i in range(0, len(line), 60)])
            #print("\n".join([line[i:i+60] for i in range(0, len(line), 60)]))
            #seq = seq + line 
            length_seq = length_seq + len(seq)     
        else:
            seq = seq + line 
            length_seq = length_seq + len(line[:-1])     

infile.close()
outfile.close()