#! /usr/bin/env python3

import numpy as np
import os
import sys
# from subprocess import *
# from Bio.PDB import *
# import openpyxl
import csv
import re
import matplotlib.pyplot as plt
import codecs
import pandas as pd
# from collections import Counter
'''26/10/17
automating conservation scores calculation for ExAC & HGMD variants'''

configfile = '/Users/mdefsss2/PycharmProjects/alignment_and_conservationScores/conservation_score_calculation1.py'
sys.path.append(os.path.dirname(os.path.expanduser(configfile)))

file = codecs.open('/Users/mdefsss2/Desktop/test_mutrot_dir/hgmd_vars.csv', "r", encoding='utf-8', errors='ignore')  # 'r' or 'rb' for binary mode and nor text only
file_reader = (csv.reader(file))
hgmd_data = list(file_reader)
file2 = codecs.open('/Users/mdefsss2/Desktop/test_mutrot_dir/training_ExAC_variants.csv', "r", encoding='utf-8', errors='ignore')
file2_reader = (csv.reader(file2))
exac_data = list(file2_reader)
file3 = codecs.open('/Users/mdefsss2/Desktop/test_mutrot_dir/hemi_gnomad_vars_20171109.csv', "r", encoding='utf-8', errors='ignore')
file3_reader = (csv.reader(file3))
gnomad_data = list(file3_reader)
file4 = pd.read_excel('/Users/mdefsss2/Desktop/test_mutrot_dir/lab_test_set_variants_analysed.xlsx')

from conservation_score_calculation1 import makeSimilarityString, alignment, symbols

## To parse the CSV file and extract res number separately in order to combine them with their conservation scores:
for seq in alignment:##  To extract the human cacna1f seq from the multi seq alignment
    # print(seq)
    # print(symbols)
    all_aa = seq[0:]
    break
# print(str(len(all_aa)))
seq_str = str(all_aa)

seq_res = []
for aa in seq_str:
    seq_res.append(aa) #  list of human residues from the alignment
# print(seq_res)
human_seq_res = []
for aa in seq_res:
    if aa != '-':
        human_seq_res.append(aa)
# print('human seq residues: ' + str(human_seq_res))
symbols_int = [int(i) for i in symbols]
# print(symbols_int)
# print(len(symbols_int))
resNumber_residue_conScore_dict = {k: v for k, v in zip(enumerate(seq_res, start=1), symbols_int) for r in range(len(seq_res)) if k[1] != '-'} ##  combining the two lists of seq residues and their corresponding conservation scores and removing the indels through '-'
# print('length of residue numbers & their conservation scores: ' + str(len(resNumber_residue_conScore_dict)))
# print('residue numbers & their conservation scores: ' + str(resNumber_residue_conScore_dict))

# modelled_res = list(key for key in resNumber_residue_conScore_dict if key[1] != '-') ##  only keeping the modelled residues but no conservation scores!
# print('length of modelled residues: ' + str(len(modelled_res)))

## testing for a mean conservation score for CTD, in the aim of comaring conservation to other regions:
# ctd_con_score = list(key for key in resNumber_residue_conScore_dict if key[0] >= 1455 and key[0] <= 1580)
# # print('ctd_con_score: ', ctd_con_score)
# exit(1)

##########################################################################
hgmd_vars = []
for row in (hgmd_data):
    hgmd_vars.append(row[1])  # column where all_hgmd_variants are
hgmd_vars.pop(0)  # delete the header in the file
# print(hgmd_vars)
print('number of all_hgmd_variants: ' + str(len(hgmd_vars)))

hgmd_res_numbers = (re.findall(r'\d+', str(hgmd_vars)))  # extracting res number
# print('residue numbers: ' + str(hgmd_res_numbers))
print('number of residues: ' + str(len(hgmd_res_numbers)))
hgmd_res_numbers_int = [int(i) for i in hgmd_res_numbers]
print('hgmd var numbers: ' + str(len(hgmd_res_numbers_int)))

# number_and_res = [k for k, v in resNumber_residue_conScore_dict.items()] ##  constructing a list of residue numbering when aligned and the residues themselves
# # print('number and residues: ' +str(number_and_res))
# res_number = [number[0] for number in number_and_res] ##  constructing a list of residue numbering as aligned in the multi alignment
# # print('aligned residue numbers: ' + str(res_number))

con_score = [v for k, v in resNumber_residue_conScore_dict.items()] ##  constructing a list of the values which are the conservation scores
# print(con_score)
resNumber_conScore_dict = {k: v for k, v in zip(enumerate(human_seq_res, start=1), con_score) for i in range(len(human_seq_res))} ##  combining the two lists of residue numbers enumerated  and their corresponding conservation scores i.e. key is made up of a tuple
# print(resNumber_conScore_dict)
seq_res_num_con_score = [k for k, v in resNumber_conScore_dict.items()] ##  constructing a list of the keys, which are made up of the human residue numbering and the residues
# print(seq_res_num_con_score)
human_seq_res_num = [num[0] for num in seq_res_num_con_score] ##  extracting the residue numbers only
# print('real res numbers: ' +str(human_seq_res_num))
human_Number_conScore_dict = {k: v for k, v in zip(human_seq_res_num, con_score) for i in range(len(human_seq_res_num))} ##  combining the two lists of the human residue numbering and their corresponding conservation scores
print('human residues & conservatin scores: ', human_Number_conScore_dict)
hgmd_human_res_con_scores = []
for number in range(len(hgmd_res_numbers_int)): ##  constructing a list to which the corresponding conservation scores are added
    hgmd_human_res_con_scores.append(human_Number_conScore_dict[hgmd_res_numbers_int[number]])
print('hgmd_human_res_con_scores: ', hgmd_human_res_con_scores)
####################################################################################################################
'''Automating con scores for ExAC variants'''

exac_vars = []
for row in (exac_data):
    exac_vars.append(row[1])  # column where all_hgmd_variants are
exac_vars.pop(0)  # delete the header in the file
# print(hgmd_vars)
print('number of ensembl_variants: ' + str(len(exac_vars)))
exac_res_numbers = (re.findall(r'\d+', str(exac_vars)))  # extracting res number
# print('residue numbers: ' + str(gnomad_numbers))
print('number of residues: ' + str(len(exac_res_numbers)))
exac_res_numbers_int = [int(i) for i in exac_res_numbers]
print('exac var numbers: ' + str(len(exac_res_numbers_int)))
exac_human_res_con_scores = []
for number in range(len(exac_res_numbers_int)): ##  constructing a list to which the corresponding conservation scores are added
    exac_human_res_con_scores.append(human_Number_conScore_dict[exac_res_numbers_int[number]])
print('exac_human_res_con_scores: ', exac_human_res_con_scores)
########################################## Automating con scores for gnomad variants:

gnomad_vars = []
for row in gnomad_data:
    gnomad_vars.append(row[1])  # column where variants are
gnomad_vars.pop(0)  # delete the header in the file
# print(gnomad_vars)
print('number of gnomad_vars: ' + str(len(gnomad_vars)))

gnomad_res_numbers = (re.findall(r'\d+', str(gnomad_vars)))  # extracting res number
# print('residue numbers: ' + str(gnomad_res_numbers))
print('number of residues: ' + str(len(gnomad_res_numbers)))
gnomad_res_numbers_int = [int(i) for i in gnomad_res_numbers]
print('gnomad var numbers: ' + str(len(gnomad_res_numbers_int)))

gnomad_human_res_con_scores = []
for number in range(len(gnomad_res_numbers_int)): ##  constructing a list to which the corresponding conservation scores are added
    gnomad_human_res_con_scores.append(human_Number_conScore_dict[gnomad_res_numbers_int[number]])
print('gnomad_human_res_con_scores: ', gnomad_human_res_con_scores)
##########################################

# lab_data = list(file4)
lab_vars = file4['Variants'].tolist()
print('number of lab_variants: ' ,len(lab_vars))
lab_res_numbers = (re.findall(r'\d+', str(lab_vars)))  # extracting res number
# print('residue numbers: ' + str(lab_numbers))
print('number of residues: ' + str(len(lab_res_numbers)))
lab_res_numbers_int = [int(i) for i in lab_res_numbers]
print('lab var numbers: ' + str(len(lab_res_numbers_int)))
lab_human_res_con_scores = []
for number in range(len(lab_res_numbers_int)): ##  constructing a list to which the corresponding conservation scores are added
    lab_human_res_con_scores.append(human_Number_conScore_dict[lab_res_numbers_int[number]])
print('lab_human_res_con_scores: ', lab_human_res_con_scores)

# n_bins = 10
myalpha = 0.3
mynormed = 1
# bin_ranges = [i for i in np.linspace(1, 9, 10)] # bin 9 & 10 are least conserved with conservation scores of 9 in the alignment (10 bins to produce a symetric histogram)
fig, ax = plt.subplots()
ax.hist(exac_human_res_con_scores, normed=mynormed, alpha=myalpha, label = "ExAC")
ax.hist(hgmd_human_res_con_scores, normed=mynormed, alpha=myalpha, label = "HGMD")
ax.hist(gnomad_human_res_con_scores, normed=mynormed, alpha=myalpha, label = "gnomAD")
ax.legend()
plt.xlabel('Conservation')
plt.ylabel('frequency')
plt.title('Distribution of Conservation Scores of Residues mutated')
fig.tight_layout()
plt.show()

