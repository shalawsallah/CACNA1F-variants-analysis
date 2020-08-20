#! /usr/bin/env python3
# 25/07/17 edited 13/12/18
#  comparing conservation scores obtained from CONSURF for gnomad & HGMD variants :
import numpy as np
import os
import sys
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd

'''changing for each gene: 
1) model = True or False 
2) variants' file 
3) score_file for each gene'''
gene = 'rpgr'
model = True
variants_file = pd.read_excel('/Users/mdefsss2/x_linked_genes/'+gene+'/'+gene+'_analysis.xlsx')
scores_file = pd.read_excel('/Users/mdefsss2/x_linked_genes/'+gene+'/consurf_scores.xlsx')

aa_format_dic = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}

vars = variants_file['variants'].tolist()
var_class = variants_file['class']
print('number of variants: ' + str(len(vars)))
res_numbers = (re.findall(r'\d+', str(vars)))
res_numbers_int = [int(i) for i in res_numbers]
print('var numbers: ' + str(len(res_numbers_int)))

consurf_score_normalized = scores_file['SCORE'].tolist()
consurf_aa_seq = scores_file[' SEQ']
consurf_colour = scores_file['COLOR'].tolist()
consurf_colour_numbers = (re.findall(r'\d+', str(consurf_colour)))
consurf_colour_score = [int(i) for i in consurf_colour_numbers]
# consurf_sequenceBased = []
# for i in consurf_colour_score:
#     consurf_sequenceBased.append(re.findall(r'\d+', str(i)))

consurf_aa_num_str = []
if model == False: #a column is missing from the generated consurf results when no pdb structure is provided
    consurf_res_num = scores_file['    3LATOM'].tolist()
    for i in consurf_res_num: #strip residue number & append
        if i != '         -':
            consurf_aa_num_str.append(re.findall(r'\d+', i))
        elif i == '         -':
            consurf_aa_num_str.append(0) #to place a digit to keep the same length as the sequence
    print('number of consurf residues: ', len(consurf_aa_num_str))
    consurf_aa_num_list = re.findall(r'\d+', str(consurf_aa_num_str)) # further strip of lists within a list
    print('number of consurf residues STILL: ', len(consurf_aa_num_list))

    # consurf_res_numbers_str = (re.findall(r'\d+', str(consurf_res_num)))
    consurf_res_numbers_int = [int(i) for i in consurf_aa_num_list] # convert to integers
    print('number of consurf residues STILL: ', len(consurf_res_numbers_int))
    # print(consurf_res_numbers_int)
    print('number of consurf residues STILL: ', len(consurf_aa_seq))

    if len(consurf_aa_seq) != len(consurf_res_numbers_int):
        print(sys.stderr)
        sys.exit('ERROR: Consurf residue names and number columns unequal when counted and stripped')

    normalized_scores = []
    colour_scores = []
    for j in range(len(vars)):
        flag= True
        for i in range(len(consurf_aa_seq)):
            # if hgmd_gene[j] == gene_of_interest:
            # print(consurf_aa_seq[i].strip())
            if aa_format_dic[vars[j].rstrip()[:3]] == consurf_aa_seq[i].strip():
                if res_numbers_int[j] == consurf_res_numbers_int[i]:
                    normalized_scores.append(consurf_score_normalized[i])
                    colour_scores.append(consurf_colour_score[i])
                    flag=False
                    # print(consurf_res_num[i], consurf_score_normalized[i])
                    # print(consurf_score_normalized[i])
                    print(consurf_colour_score[i])
        if flag==True:
            print(' ')
    # print('hgmd scores: ', hgmd_scores)

else: #model & not pdb structure
    consurf_res_numbers_int = scores_file[' POS'].tolist()
    print('number of consurf residues STILL: ', len(consurf_res_numbers_int))
    print(consurf_res_numbers_int)
    print('number of consurf residues STILL: ', len(consurf_aa_seq))

    if len(consurf_aa_seq) != len(consurf_res_numbers_int):
        print(sys.stderr)
        sys.exit('ERROR: Consurf residue names and number columns unequal when counted and stripped')

    normalized_scores = []
    colour_scores = []
    for j in range(len(vars)):
        flag= True
        for i in range(len(consurf_aa_seq)):
            # if hgmd_gene[j] == gene_of_interest:
            # print(consurf_aa_seq[i].strip())
            if aa_format_dic[vars[j].rstrip()[:3]] == consurf_aa_seq[i].strip():
                if res_numbers_int[j] == consurf_res_numbers_int[i]:
                    normalized_scores.append(consurf_score_normalized[i])
                    colour_scores.append(consurf_colour_score[i])
                    flag=False
                    # print(consurf_res_numbers_int[i], consurf_aa_seq[i].rstrip(), consurf_score_normalized[i])
                    # print(consurf_score_normalized[i])
                    print(consurf_colour_score[i])
        if flag==True:
            print(' ')
    # print('hgmd scores: ', hgmd_scores)

    '''getting rid of the * present with some of the scores (calculated from less than 6 non-gapped sequences) in consurf 'COLOR' '''
    # hgmdScores = (re.findall(r'\d+', str(hgmd_colour_scores)))
    # gnomadScores = (re.findall(r'\d+', str(gnomad_colour_scores)))
    # hgmdScores_int = [int(i) for i in hgmdScores]
    # gnomadScores_int = [int(i) for i in gnomadScores]

# pathogenic_consurf_scores = []
# nonpathogenic_consurf_scores = []
# for i in range(len(vars)):
#     if var_class[i] == 'pathogenic':
#         pathogenic_consurf_scores.append(consurf_colour_score[i])
#     elif var_class[i] == 'benign':
#         nonpathogenic_consurf_scores.append(consurf_colour_score[i])
#
# mynormed = 1
# myalpha = 0.5
# fig, ax = plt.subplots()
# ax.hist(nonpathogenic_consurf_scores, normed=mynormed, histtype='step', stacked=True, fill=False, alpha=myalpha,
#         label="Dataset N")
# ax.hist(pathogenic_consurf_scores, normed=mynormed, histtype='step', stacked=True, fill=False, alpha=myalpha,
#         label="Dataset P")
# # ax.hist(gnomad_scores, bin_ranges, normed=1, histtype='step', stacked=True, fill=False, alpha=myalpha, label = "gnomAD")
# ax.legend()
# plt.xlabel('Conservation')
# plt.ylabel('Frequency')
# # plt.title('Consurf Conservation Scores (normalized) of Residues mutated')
# fig.tight_layout()
# plt.show()

# fig, ax = plt.subplots()
# ax.hist(gnomadScores_int, normed=mynormed, histtype='step', stacked=True, fill=False, alpha=myalpha, label="gnomAD")
# ax.hist(hgmdScores_int, normed=mynormed, histtype='step', stacked=True, fill=False, alpha=myalpha, label="HGMD")
# # ax.hist(gnomad_scores, bin_ranges, normed=1, histtype='step', stacked=True, fill=False, alpha=myalpha, label = "gnomAD")
# ax.legend()
# plt.xlabel('Conservation')
# plt.ylabel('frequency')
# plt.title('Consurf Conservation Scores of Residues mutated')
# fig.tight_layout()
# plt.show()


