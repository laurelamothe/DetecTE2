<<<<<<< HEAD
##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                                          TE trimming :                                  #############
############                      a script to trimm the sequence shorter than 500 bp                 #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
parser = argparse.ArgumentParser(prog="TE_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--outfile", "-o", required=True, type=str, action="store", help="Path to the output file.")
args = parser.parse_args()

infile = open(args.fasta, "r")
outfile = open(args.outfile, "w")
title=""          
seq=""  
length_seq = 0        
for line in infile:
    if line[0] == '>':
        if length_seq >= 500:
            outfile.write(title)
            outfile.write(seq)
        title = line
        seq=""
        length_seq = 0
    else:
        seq = seq + line 
        length_seq = length_seq + len(line[:-1])     

infile.close()
=======
##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                                          TE trimming :                                  #############
############                      a script to trimm the sequence shorter than 500 bp                 #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
parser = argparse.ArgumentParser(prog="TE_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--outfile", "-o", required=True, type=str, action="store", help="Path to the output file.")
args = parser.parse_args()

infile = open(args.fasta, "r")
outfile = open(args.outfile, "w")
title=""          
seq=""  
length_seq = 0        
for line in infile:
    if line[0] == '>':
        if length_seq >= 500:
            outfile.write(title)
            outfile.write(seq)
        title = line
        seq=""
        length_seq = 0
    else:
        seq = seq + line 
        length_seq = length_seq + len(line[:-1])     

infile.close()
>>>>>>> 240ff64487a463196faacf2189db5ea1f6c2dcf5
outfile.close()