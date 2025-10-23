#!/usr/bin/python3
import argparse
import os.path

parser = argparse.ArgumentParser(description='Create file with the consensus sequence keeping insertions')
parser.add_argument('-f', '--fasta', dest='fastaFile', type=str, help='input fasta file of a multiple alignment', required=True)
parser.add_argument('-o','--output', dest='output_file', type=str, default="consensus_withInserts.fa", help="output file")
args = parser.parse_args()

dic={}

#Create dictionary
fil=open(args.fastaFile, "r")
for line in fil.readlines() :
    if str(line)[0] == '>':
        fid=str(line.splitlines()[0])
        dic.setdefault(fid, [])
    else:
        for car in str(line):
            if car!='\n' :
                dic[fid].append(car)
tail_seq=len(dic[fid])
fil.close()

#Get consensus
i=0
consensus=""
while i < tail_seq :
    tokeep={}
    for el in dic:
        if dic[el][i] in tokeep.keys():
            tokeep[dic[el][i]]=tokeep[dic[el][i]]+1
        else:
            tokeep[dic[el][i]]=1
    if "-" in tokeep.keys():
        del tokeep["-"]
    if len(tokeep) > 1 :
        sorted_tokeep= sorted(tokeep.items(), key=lambda x: x[1], reverse=True)
        toprint= sorted_tokeep[0][0]
        if sorted_tokeep[0][1] == sorted_tokeep[1][1] :
            if sorted_tokeep[0][0] == '-' :
                toprint= sorted_tokeep[1][0]
#		print(sorted_tokeep)
    else :
        toprint=list(tokeep.keys())[0]
    consensus=consensus+str(toprint)
    i=i+1

seqname=os.path.basename(args.fastaFile).split(".fa")
seqN=seqname[0]
f = open(args.output_file, "w")
f.write(">" + seqN + "\n")
f.write(consensus + "\n")
f.close()
