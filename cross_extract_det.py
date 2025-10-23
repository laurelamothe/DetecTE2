##################################################################################################################
##################################################################################################################
##################################################################################################################
############                                                                                         #############
############                          det_cross_extract :  a script to select                        #############
############               the transcripts with a crossmatch between blastx and tblastn              #############
############                                                                                         #############
##################################################################################################################
##################################################################################################################
##################################################################################################################

import argparse
import pandas as pd
import numpy as np
from itertools import repeat

parser = argparse.ArgumentParser(prog="transcript_size_filter.py", description="write a fasta with only the sequences longer than $SIZE nucl")
parser.add_argument("--fasta", "-f", required=True, type=str, action="store", help="Nucl fasta file")
parser.add_argument("--blastx", "-b", required=True, type=str, action="store", help="Path to blastx output file")
parser.add_argument("--tblastn", "-t", required=True, type=str, action="store", help='Path to tblastn output file')
parser.add_argument("--outdir", "-o", required=True, type=str, action="store", help="Path to the output directory. The directory will be created during execution. If not specified, the directory will be created in the script location.")
parser.add_argument("--prefix", "-p", required=True, type=str, action="store", help='Prefix for the output files. If not specified it will take the letters before "." in fasta file name')
args = parser.parse_args()

#V#######################################################################################
#################                       CROSSMATCH                      #################
#########################################################################################

# Checks that blastx and tblastn show consistent results
# output : 
#       - list_blastx / list_tblastn = list of blast output lines + corresonding TE order and superfamilly
#       - list_cross_O / list_cross_SF = list of blast output lines meeting crossmatch requirement at Order or Superfamilly level
#       - dico_blastx_O / dico_tblastn_O / dico_blastx_SF / dico_tblastn_SF = dico of sequence meetin crossmatch requirement. key is sequence ID transcriptome, value is the order or super familly 

list_blastx = []
infile = open(args.blastx, "r")
for line in infile:
    item = line[:-1].split("\t")
    ID = item[0]
    annot = item[1].split("#")
    annot2 = annot[1].replace("-", "/").split("/")
    new_line = [ID] + [annot2[0], annot2[2]]
    list_blastx.append(new_line)

list_tblastn = []
infile = open(args.tblastn, "r")
for line in infile:
    item = line[:-1].split("\t")
    ID = item[1]
    annot = item[0].split("#")
    annot2 = annot[1].replace("-", "/").split("/")
    new_line = [ID] + [annot2[0], annot2[2]]
    list_tblastn.append(new_line)
  
dico_blastx_SF = {}  # {key = ID : value = SUPERFAMILLY}
for hit in list_blastx:
    dico_blastx_SF[hit[0]] = hit[-1]

dico_tblastn_SF = {}  # {key = ID : value = list_of_SUPERFAMILLIES}
for hit in list_tblastn:
    if hit[0] not in dico_tblastn_SF:
        dico_tblastn_SF[hit[0]] = [hit[-1]]
    else:
        dico_tblastn_SF[hit[0]].append(hit[-1])


#### CROSSMATCH ON SUPERFAMILLY ####
list_cross_SF = list(list_blastx)
for element in list_cross_SF:
    id = element[0]
    if id in dico_tblastn_SF:
        if dico_blastx_SF[id] not in dico_tblastn_SF[id]:
            list_cross_SF.remove(element)
    else:
        list_cross_SF.remove(element)
        
list_transcrit_cross=[]
for element in list_cross_SF:
    list_transcrit_cross.append(element[0])


#### WRITE TBLASTN OUTPUT WITH CROSSMATCHED TRANSCRIPTS ####
infile = open(args.tblastn, "r")
outfile = open("{}{}_TBLASTN_cross.out".format(args.outdir, args.prefix), "w")
for line in infile : 
    tmp = line.split("\t")
    if tmp[1] in list_transcrit_cross :
        outfile.write(line)
outfile.close()



#########################################################################################
#################                       EXTRACTION                      #################
#########################################################################################

# Extract the sequence of transposable elements in the crossmatched transcripts from tblastn output and annotate them
#   1) merge the hits of same superfamilly
#       - if two hits are separated by more than 250 bp, there are considered to arise from two distinct elements
#       - all e-values are stored, as is the superfamilly, the transcript name, the direction and limits of each hit
#       - output : tab1 ['transcript', 'annot', 'dir', 'positions', 'limits', 'e-values']
#   2) merge the putative TE of different superfamillies
#       - if two putative TE are not overlapping, they are considered to arise from two distinct elements
#       - output : tab3 ['transcript', 'annot', 'dir', 'positions', 'limits', 'e-values']
#           . with one line per confirmed TE, with list of possible anotation, all e-values, all positions of the putative TE and extern limits of the TE
#   3) find the gaps in the confirmed TEs
#       - add column 'gaps' and 'gap_range' to tab3
#       - search for the intersections in the gaps of each putative TE
#   4) Annotate the confirmed TEs
#       - create tab4 ['annot', 'e-value'] with every hit of the the confirmed TE
#       - sort the hits by e-values
#       - in the 10 first hits (or less), if 80% have the same annotation, the TE is annotated, else it id anotated (ODER/Unknown of Unknown/Unknown)


#### MERGE HITS PER SUPERFAMILLIES ####
# Edit the tblastn output to add direction and format indexes
all_cross = pd.read_table("{}{}_TBLASTN_cross.out".format(args.outdir, args.prefix), sep="\t",
                          names = ['query_seq_id', 'db_seq_id', 'percent_ident', 'length', 'mismatch', 'gap_open', 'query_start', 'query_end', 'db_start', 'db_end', 'e_value', 'bitscore'] )
SF_list = []
for line in all_cross['query_seq_id'] :
    SF_list.append(line.split('#')[1].split("-")[0])
all_cross['query_seq_id'] = SF_list # replace db seq name by anotation
list_transcrit_cross= np.unique(all_cross['db_seq_id'])
list_dir = list(repeat("+", len(all_cross.index)))
all_cross.db_start = np.array(all_cross.db_start)-1
all_cross.db_end = np.array(all_cross.db_end)-1
for line in all_cross.index : 
    if all_cross['db_start'][line] > all_cross['db_end'][line] : # update the direction and invert start and end indexes if  TE in compl strand
        list_dir[line] = '-' 
        tmp = all_cross['db_start'][line]
        all_cross._set_value(line, 'db_start', all_cross['db_end'][line])
        all_cross._set_value(line, 'db_end', tmp)
all_cross.insert(10, "dir", list_dir)

# write tab1 with one ligne per transcript and superfamilly 
tab1  = pd.DataFrame({'transcript' : [], 'annot' : [], 'dir' : [], 'positions' : [], 'limits' : [], 'e-values' : []})
for transcript in np.unique(all_cross['db_seq_id']) : 
    current_transcrip_table = all_cross[all_cross['db_seq_id'] == transcript] # select all hits 
    for SF in np.unique(list(current_transcrip_table['query_seq_id'])) :
        for direction in ['+', '-'] : 
            current_tblastn_table = current_transcrip_table[current_transcrip_table['query_seq_id'] == SF] # select one superfamilly
            current_tblastn_table = current_tblastn_table[current_tblastn_table['dir'] == direction] # select one direction
            e_values = []
            if len(current_tblastn_table) > 0 : 
                current_tblastn_table = current_tblastn_table.sort_values(by=['db_start','db_end'], ascending=[True, True]) # sort hits by position
                current_lims = [[current_tblastn_table['db_start'][current_tblastn_table.index[0]], current_tblastn_table['db_end'][current_tblastn_table.index[0]]]] # store position of first hit
                if len(current_tblastn_table) > 1 :
                    for hit in range(0, len(current_tblastn_table)-1) :
                        e_values.append(current_tblastn_table['e_value'][current_tblastn_table.index[hit]]) # store e-values
                        next_lims = [current_tblastn_table['db_start'][current_tblastn_table.index[hit + 1]], current_tblastn_table['db_end'][current_tblastn_table.index[hit + 1]]] # store positions next hit
                        gap  = next_lims[0] - current_lims[-1][1]
                        if gap < 1 : # if overlapping, update current limits
                            current_lims[-1] = [min([current_lims[-1][0], next_lims[0]]), max([current_lims[-1][1], next_lims[1]])] 
                        elif gap < 250 : # if small gap, add next to current
                            current_lims.append(next_lims) 
                        else : # if big gap, write current in tab1 (new putative TE!) and start new TE
                            tab1 = tab1._append({'transcript' : current_tblastn_table['db_seq_id'][current_tblastn_table.index[hit]], 'annot' : current_tblastn_table['query_seq_id'][current_tblastn_table.index[hit]], 'dir' : direction, 'positions' : current_lims, 'limits' : [current_lims[0][0], current_lims[-1][1]], 'e-values' : e_values}, ignore_index = True) 
                            current_lims = [next_lims] 
                            e_values = []
                e_values.append(current_tblastn_table['e_value'][current_tblastn_table.index[-1]]) 
                tab1 = tab1._append({'transcript' : current_tblastn_table['db_seq_id'][current_tblastn_table.index[-1]], 'annot' : current_tblastn_table['query_seq_id'][current_tblastn_table.index[-1]], 'dir' : direction, 'positions' : current_lims, 'limits' : [current_lims[0][0], current_lims[-1][1]], 'e-values' : e_values}, ignore_index = True) # write final putative TE in tab1 

                    
#### MERGE HITS ON SAME TRANSCRIPT ####
tab3  = pd.DataFrame({'transcript' : [], 'annot' : [], 'dir' : [], 'positions' : [], 'limits' : [], 'e-values' : []})
for transcript in np.unique(tab1['transcript']) : 
    for direction in ['+', '-'] : 
        tab2 = tab1.loc[tab1['transcript'] == transcript] # select putative TE on 1 transcript
        tab2 = tab2.loc[tab2['dir'] == direction] # select one direction
        if len(tab2) > 0 :
            tab2 = tab2.sort_values(by = ['limits'])
            current_lim = tab2['limits'][tab2.index[0]] # store current putative TE limits
            current_pos = [tab2['positions'][tab2.index[0]]] # store all hits
            current_SF = [tab2['annot'][tab2.index[0]]] # store anotation
            current_evalue = [tab2['e-values'][tab2.index[0]]] # store e-value
            if len(tab2) > 1 :
                next_lim = []
                for line in range(0, len(tab2) - 1) : 
                    next_lim = tab2['limits'][tab2.index[line + 1]] # store next putative TE lims
                    if current_lim[1] >= next_lim[0] : # if overlapping, update putative TE infos
                        current_lim = [min([current_lim[0], next_lim[0]]), max([current_lim[1], next_lim[1]])]
                        current_pos.append(tab2['positions'][tab2.index[line + 1]])
                        current_SF.append(tab2['annot'][tab2.index[line + 1]])
                        current_evalue.append(tab2['e-values'][tab2.index[line + 1]])
                    else : # if gap, write new confirmed TE in tab3 and start new TE
                        tab3 = tab3._append({'transcript' : tab2['transcript'][tab2.index[line]], 'annot' : current_SF, 'dir' : direction, 'positions' : current_pos, 'limits' : current_lim, 'e-values' : current_evalue}, ignore_index = True)
                        current_lim = tab2['limits'][tab2.index[line + 1]]
                        current_pos = [tab2['positions'][tab2.index[line + 1]]]
                        current_SF = [tab2['annot'][tab2.index[line + 1]]]
                        current_evalue = [tab2['e-values'][tab2.index[line + 1]]]
            tab3 = tab3._append({'transcript' : tab2['transcript'][tab2.index[0]], 'annot' : current_SF, 'dir' : direction, 'positions' : current_pos, 'limits' : current_lim, 'e-values' : current_evalue}, ignore_index = True) # write final confirmed TE
# tab3 now includes all confirmed TE


#### FIND GAPS IN CONFIRMED TE ####
tab3.insert(4, 'gaps', list(repeat([], len(tab3))))
for index in tab3.index : 
    lims = tab3['limits'][index]
    gap_list = []
    for SF in range(0, len(tab3['annot'][index])) :
        range_list = []
        if  lims[0] < tab3['positions'][index][SF][0][0] : # add gap between TE start and first putative TE piece start in the current superfamilly
            range_list .append(list(range(lims[0], tab3['positions'][index][SF][0][0]+1)))
        if tab3['positions'][index][SF][-1][1] < lims[1] : # idem on TE end
            range_list.append(list(range(tab3['positions'][index][SF][-1][1], lims[1]+1)))
        if len(tab3['positions'][index][SF]) > 1 : 
            for element in range(0,len(tab3['positions'][index][SF])-1) : # store each gap beteween putative TE pieces
                range_list.append(list(range(tab3['positions'][index][SF][element][1], tab3['positions'][index][SF][element+1][0]+1)))
        if len(range_list) > 0 :
            range_list = np.unique(np.concatenate(range_list))
        gap_list.append(set(range_list))
        final_gaps =np.unique(gap_list[0].intersection(*gap_list)) # intersection of gaps of all putative TE = confirmed TE gaps
        final_gaps.sort()
    tab3['gaps'][index] = list(gap_list[0].intersection(*gap_list)) # update tab3 

tab3.insert(6, 'gap_range', list(repeat([], len(tab3))))
for index in range(0, len(tab3.gaps)) :  # add gaps lims (just output sequence titles)             
    gap_list = []         
    a = tab3.gaps[index] 
    if a != [] : 
        lim = [a[0], a[0]]
        for i in range(0,len(a[:-1])) : 
            if a[i+1] - a[i] > 1 : 
                lim[1] = a[i]
                gap_list.append(lim)
                lim = [a[i+1], a[i+1]]
        lim[1] = a[-1]
        gap_list.append(lim)
    tab3.gap_range[index] = gap_list

   
#### ANNOTATE CONFIRMED TE ####
tab3.insert(2, 'final_annot', list(repeat([], len(tab3))))
for index in tab3.index :
    tab4 = pd.DataFrame({'annot' : [], 'e-value' : []}) 
    for sf in range(0, len(tab3['annot'][index])) : # write tab4 with every hit's annotation and e-value
        SF = tab3['annot'][index][sf]
        for e_value in tab3['e-values'][index][sf] : 
 
            tab4 = tab4._append({'annot' : SF, 'e-value' : e_value}, ignore_index = True)
    tab4 = tab4.sort_values(by = 'e-value')
    tab4 = tab4.reset_index()[['annot', 'e-value']]
    sub_tab4 = tab4.loc[tab4.index < 10] # select 10 best hits
    best_oders = np.unique([sub_tab4['annot'][i].split("/")[0] for i in range(0, len(sub_tab4['annot']))])
    if len(best_oders) > 1 : # if several orders in bests hist -> Unknown
        final_annot = "InterOrderUnknown-" + "-".join(list([i.split("/")[0] for i in np.sort(np.unique(sub_tab4['annot']))])) + "/Unknown"
    else : 
        final_annot = best_oders[0] + '/Unknown'
        for annot in np.unique(sub_tab4['annot']) : # if 80% of best hits are consistent -> annotate
            if len(sub_tab4.loc[sub_tab4['annot'] == annot]) >= 0.8 * len(sub_tab4) : 
                final_annot = annot
        if final_annot == best_oders[0] + '/Unknown' : 
            best_group = np.unique([i.split("/")[1] for i in np.unique(sub_tab4['annot'])])
            if  len(best_group) > 1 : 
                final_annot = best_oders[0] + '/InterGroupUnknown-' + "-".join(list([i.split("/")[1] for i in np.sort(np.unique(sub_tab4['annot']))]))
            else : 
                final_annot = best_oders[0] + '/' + best_group[0] + '/IntraGroupUnknown-' + "-".join(list([i.split("/")[2] for i in np.sort(np.unique(sub_tab4['annot']))]))
    tab3['final_annot'][index] = final_annot # update tab3


#VI######################################################################################
#################                      OUTPUT FILES                     #################
#########################################################################################

# creates the final output files
# outputs : 
#       - dico_cross_SF = dictionary of transcripts confirmed by crossmatch sorted by superfamilly
#       - [prefix]_confirmedTE_<ORDER>-<SUPERFAMILLY>.fasta = .fasta file with confirmed TE sorted by superfamilly
#       - [prefix]_confirmedTE.fasta = .fasta file with all confirmed TE


#### FASTA WITH CROSSMATCHED TRANSCRIPTS ####
infile =  open(args.fasta, "r")
outfile = open("{}{}_cross_transcripts.fasta".format(args.outdir, args.prefix), "w")
for line in infile :
    if line[0] == '>' :  
        if line[1:-1] in list_transcrit_cross :
            flag = True
        else : 
            flag = False
    if flag : 
        outfile.write(line)
outfile.close()

# crossmatched transcripts with seq in one line
infile = open("{}{}_cross_transcripts.fasta".format(args.outdir, args.prefix), "r")
outfile = open("{}/tmp/{}_cross_transcripts_oneline.fasta".format(args.outdir, args.prefix), "w")
found=False
c=0
for line in infile: 
    if line[0] == '>': 
        if line[1:-1] in list(tab3['transcript']):
            found = True
            if c == 0: 
                line = line + "\n"
                c+=1
            else :
                line = "\n" + line + "\n"
        else : 
            found = False               
    if found:
        outfile.write(line[:-1])
outfile.close()
infile.close()

#### FASTA WITH ALL CONFIRMED TE ####
# correspondence key for reverse complement
comp = {
    'A' : 'T',
    'T' : 'A',
    'G' : 'C',
    'C' : 'G',
    'N' : 'N'}

infile = open("{}/tmp/{}_cross_transcripts_oneline.fasta".format(args.outdir, args.prefix), "r")
outfile = open("{}{}_confirmedTE.fasta".format(args.outdir, args.prefix), "w")
for line in infile : 
    if line[0] == '>' : # store transcript
        current_transcript_table = tab3.loc[tab3['transcript'] == str(line[1:-1])]
        current_transcript_table =  current_transcript_table.reset_index().iloc[:, 1:]
    else : 
        for TE in range(0, len(current_transcript_table)) : 
            outfile.write('>' + current_transcript_table['transcript'][TE] + '_' + str(TE) + '#' + current_transcript_table['final_annot'][TE] + ' : ' + str(current_transcript_table['limits'][TE]) + ', direction = ' + current_transcript_table['dir'][TE] + ', gaps = ' + str(current_transcript_table['gap_range'][TE]) + '\n') # write sequence title (transcript_n°#ORDER/Superfamilly : [limits], direction = '', gaps = [])
            seq= ''
            for nucl in range(0, len(line)-1) : # extract confirme TE sequence in '+' strand
                if current_transcript_table['limits'][TE][0] <= nucl <= current_transcript_table['limits'][TE][1] and nucl not in current_transcript_table['gaps'][TE] : 
                    seq = "".join([seq, line[nucl]])
            if current_transcript_table.dir[TE] == '-' : # reverse complement if needed
                rev_seq = "".join(reversed(seq))
                comp_seq = ''
                for nucl in rev_seq : 
                    comp_seq = "".join([comp_seq, comp[nucl]])
                seq = comp_seq
            outfile.write("\n".join([seq[i:i+60] for i in range(0, len(seq), 60)]) + "\n") # write seq in 60 nucl increments
infile.close()
outfile.close()          
