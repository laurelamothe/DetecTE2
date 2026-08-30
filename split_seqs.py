<<<<<<< HEAD
##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                            Split Seqs : a script cut a fasta                            #############
############                    in several file with a given number of sequences                     #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
parser = argparse.ArgumentParser(prog="transcript_sizr_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help='Prefix for the output files. If not specified it will take the letters before "." in fasta file name')
parser.add_argument("--seq_num", "-n", required=True, type=int, action="store", help='Number of sequence in each output fasta file')
args = parser.parse_args()


infile = open(args.fasta, "r")
suf = 1
outfile = open("{}{}-{}.fasta".format(args.outdir, args.prefix, suf), "w")
c = 0
for line in infile : 
    if line[0] == '>' : 
        c+=1
    if c >= args.seq_num :
        c=0
        outfile.close()
        suf+=1
        outfile=open("{}{}-{}.fasta".format(args.outdir, args.prefix, suf), "w") 
    outfile.write(line)
outfile.close()
infile.close()
=======
##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                            Split Seqs : a script cut a fasta                            #############
############                    in several file with a given number of sequences                     #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
parser = argparse.ArgumentParser(prog="transcript_sizr_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")

parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help='Prefix for the output files. If not specified it will take the letters before "." in fasta file name')
parser.add_argument("--seq_num", "-n", required=True, type=int, action="store", help='Number of sequence in each output fasta file')
args = parser.parse_args()


infile = open(args.fasta, "r")
suf = 1
outfile = open("{}{}-{}.fasta".format(args.outdir, args.prefix, suf), "w")
c = 0
for line in infile : 
    if line[0] == '>' : 
        c+=1
    if c >= args.seq_num :
        c=0
        outfile.close()
        suf+=1
        outfile=open("{}{}-{}.fasta".format(args.outdir, args.prefix, suf), "w") 
    outfile.write(line)
outfile.close()
infile.close()
>>>>>>> 240ff64487a463196faacf2189db5ea1f6c2dcf5
