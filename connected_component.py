##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                         connected_component :  a script to merge                        #############
############                the identical TE from different assemblies of a transcriptome            #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
import networkx as nx
import pandas as pd
import numpy as np
import itertools as it

parser = argparse.ArgumentParser(prog="transcript_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")
parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--blast", "-b", required=True, type=str, action="store", help="Path to blastx output file")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help='Prefix for the output files. If not specified it will take the letters before "." in fasta file name')
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
args = parser.parse_args()

#########################################################################################
#################                 TREATMENT OF UNKNOWN                  #################
#########################################################################################
# removes the unknowns that blasted with non-unknown sequences to avoid 'bridges' between famillies
df = pd.read_csv(args.blast, index_col=False, skiprows=0, sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
df2 = df[df.pident >= 90]
df2 = df2[df2.length > 100].reset_index(drop=True, inplace=False)
seq_to_drop = [df2.qseqid[i] for i in df2.index if "Unknown" in df2.qseqid[i]  and "Unknown" not in df2.sseqid[i]]
seq_to_drop = seq_to_drop + [df2.sseqid[i] for i in df2.index if "Unknown" in df2.sseqid[i]  and "Unknown" not in df2.qseqid[i]]
seq_to_drop = np.unique(seq_to_drop)
lines_to_drop = [i for i in df2.index if df2.qseqid[i] in seq_to_drop  or df2.sseqid[i] in seq_to_drop]
df2 = df2.drop(lines_to_drop , axis=0)

#########################################################################################
#################          TREATMENT OF MONOSEQ FAMILLIES               #################
#########################################################################################
# for the seq that, with the current filter, remains alone, adds the hits with 95% id and more than 100 bp
# monoseq = [ i for i in np.unique(df2.qseqid) if len(df2[df2.sseqid == i]) == 1]
# a=[]
# for i in df.index : 
#     if (df.sseqid[i] in monoseq or df.qseqid[i] in monoseq) and int(df.length[i]) > 100 and int(df.length[i]) <= 200 and df.pident[i] >= 90:
#         df2 = df2._append(df.iloc[i,:])
#         if df.sseqid[i] in monoseq :
#               a.append(df.sseqid[i])  
#         if df.qseqid[i] in monoseq :
#             a.append(df.qseqid[i])
# a = np.unique(a)

#########################################################################################
#################                       CLUSTERING                      #################
#########################################################################################
# groups the sequences with 99% of identity 2 by 2

#### GRAPH ####
G=nx.from_pandas_edgelist(df2, source='qseqid', target='sseqid', edge_attr=None, create_using=None, edge_key=None)
cc=sorted(nx.connected_components(G), key=len)


#### NEW ANNOTATION ####
anotated_clusters = {}
for component in cc : 
    anotations = np.unique([i.split('#')[1] for i in component])
    orders = np.unique([i.split('/')[0] for i in anotations])
    groups = np.unique([i.split('/')[1] for i in anotations])
    if len(orders) > 1 : 
        new_anot = 'Unknown/Unknown'
    elif len(groups) > 1 :         
        mix = list(filter(lambda n: n != 'Unknown' ,list(it.chain.from_iterable([i.split('/')[1].split("-") for i in component])))) # merge all annotations and remove 'Unknown' anotations
        new_anot = orders[0] + '/InterGroupMix-' + "-".join(np.unique(mix))
    elif len(anotations) > 1 :
        mix = list(filter(lambda n: n != 'Unknown' ,list(it.chain.from_iterable([i.split('/')[2].split("-") for i in component])))) # merge all annotations and remove 'Unknown' anotations
        new_anot = orders[0] + '/' + groups[0] + '/IntraGroupMix-' + "-".join(np.unique(mix))
        # for anot in np.unique(mix) : 
        #     if len(list(filter(lambda n: n == anot, mix))) >= 0.8 * len(mix) : 
        #         new_anot = orders[0] + '/' + anot   
    else : 
        new_anot = anotations[0]
    existing=[]
    for anot in anotated_clusters : 
        if anot.split('_')[0] == new_anot : 
            existing.append(int(anot.split('_')[1]))
    if len(existing) != False : 
        num_clus = max(existing) + 1
    else : 
        num_clus = 1
    anotated_clusters[new_anot.split(" ")[0] + '_' + str(num_clus)] = list(component)

#########################################################################################
#################                         OUTPUTS                       #################
#########################################################################################
# - fasta by cluster
# - table with familly counts

superfamilies = list(pd.read_csv("superfamily_list.txt", index_col=False, skiprows=0, names=['Annot']).Annot)


#### FASTA BY FAMILLY ####
seq_table = pd.DataFrame({'seq' : [], 'end_title' : []})
infile = open(args.fasta, 'r')
title = ''
end_title = ''
seq = ''
for line in infile : 
    if line[0] == '>' :
        seq_table = pd.concat([seq_table, pd.DataFrame({'seq' : [seq] , 'end_title' : [end_title]}, index=[title])])
        title = line[1:].split(' ')[0]
        end_title = ' '.join(line[1:].split(' ')[1:])
        seq = ''
    else : 
        seq += line
seq_table = seq_table.iloc[seq_table.index != '']   
infile.close()

for cluster in anotated_clusters : 
    component_seqs = seq_table.filter(items=anotated_clusters[cluster], axis=0)
    outfile = open(args.outdir + '-'.join(cluster.split('/')) + '.fasta', 'w')
    for seq in component_seqs.index : 
        outfile.write('>' + seq + '_' + cluster.split('_')[1] + ' ' + component_seqs.end_title[seq] + component_seqs.seq[seq])
    outfile.close()


#### COUNT TABLE ####
famillies = []

for i in anotated_clusters : 
    code = i.split("_")[0]
    # if "Intra" in i : 
    #     famillies.append("/".join(code.split("/")[0:2]) + "/Unknown")
    # else : 
    famillies.append("/".join(code.split("/")[0:2]) + "/Total")
    famillies.append("/".join(code.split("/")[0:1]) + "/Total/-/")
    famillies.append(code)

count_table = pd.DataFrame({"Order" : [], "Group" : [], "Superfamilly" : [], args.prefix : []})
for sf in superfamilies : 
    # fam_count = len(list(filter(lambda fam: fam == sf, famillies)))
    fam_count = famillies.count(sf)
    count_table = count_table._append({"Order" : sf.split("/")[0], "Group" : sf.split("/")[1], "Superfamilly" : sf.split("/")[2], args.prefix : fam_count}, ignore_index=True)
count_table.T.to_csv("{}count_table.csv".format(args.outdir), header = False)
