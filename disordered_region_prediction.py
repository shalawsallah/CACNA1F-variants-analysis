#! /usr/bin/env python3
import numpy as np
import os
import sys
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd

'''predicting the disordered regions from IUPRED2A @ https://iupred2a.elte.hu
 (will only include one prediction, ie uipred and not anchore, as both correlate highly) 

changing for each gene: 
1) iupred predictions or anchor 
2) variants' file 
3) predictions file for each gene'''
gene = 'rpgr'
iupred_prediction = True
anchor_prediction = False
variants_file = pd.read_excel('/Users/mdefsss2/x_linked_genes/'+gene+'/'+gene+'_analysis.xlsx')
scores_file = pd.read_excel('/Users/mdefsss2/x_linked_genes/'+gene+'/iupred2a_predictions.xlsx')

aa_format_dic = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}

vars = variants_file['variants'].tolist()
print('number of variants: ', len(vars))
res_numbers = (re.findall(r'\d+', str(vars)))
res_numbers_int = [int(i) for i in res_numbers]
print('var numbers: ' + str(len(res_numbers_int)))

predictions1 = scores_file['IUPRED SCORE'].tolist()
predictions2 = scores_file['ANCHOR SCORE'].tolist()
aa_seq = scores_file['AMINO ACID'].tolist()
aa_position = scores_file['# POS'].tolist()

predictions = []
if iupred_prediction == True:
    for k in predictions1: #strip residue number & append
        predictions.append(k)
    print('number of protein residues: ', len(predictions))
    # aa_num_list = re.findall(r'\d+', str(aa_num)) # further strip of lists within a list
    # print('number of protein residues STILL: ', len(aa_num_list))
    # aa_num_int = [int(i) for i in aa_num_list] # convert to integers
    # print('number of protein residues STILL: ', len(aa_num_int))
    # print(aa_num_int)
    # print('number of protein residues STILL: ', len(aa_seq))

    if len(aa_seq) != len(predictions):
        print(sys.stderr)
        sys.exit('ERROR: Protein residue names and number columns unequal when counted and stripped')

    scores = []
    for j in range(len(vars)):
        flag = True
        for i in range(len(aa_seq)):
            # if hgmd_gene[j] == gene_of_interest:
            # print(consurf_aa_seq[i].strip())
            if aa_format_dic[vars[j].rstrip()[:3]] == aa_seq[i].strip():
                if res_numbers_int[j] == aa_position[i]:
                    scores.append(predictions1[i])
                    flag = False
                    # print(consurf_res_num[i], consurf_score_normalized[i])
                    # print(consurf_score_normalized[i])
                    print(predictions1[i])
            elif aa_format_dic[vars[j].rstrip()[:3]] != aa_seq[i].strip():
                if res_numbers_int[j] == aa_position[i]:
                    flag = False
                    print(' ')
        if flag == True:
            print(vars[j])
elif anchor_prediction == True:
    for k in predictions2: #strip residue number & append
        predictions.append(k)
    print('number of protein residues: ', len(predictions))
    # aa_num_list = re.findall(r'\d+', str(predictions)) # further strip of lists within a list
    # print('number of protein residues STILL: ', len(aa_num_list))

    # consurf_res_numbers_str = (re.findall(r'\d+', str(consurf_res_num)))
    # aa_num_int = [int(i) for i in aa_num_list] # convert to integers
    # print('number of protein residues STILL: ', len(aa_num_int))
    # print(aa_num_int)
    # print('number of protein residues STILL: ', len(aa_seq))

    if len(aa_seq) != len(predictions):
        print(sys.stderr)
        sys.exit('ERROR: Protein residue names and number columns unequal when counted and stripped')

    scores = []
    for j in range(len(vars)):
        flag = True
        for i in range(len(aa_seq)):
            # if hgmd_gene[j] == gene_of_interest:
            # print(consurf_aa_seq[i].strip())
            if aa_format_dic[vars[j].rstrip()[:3]] == aa_seq[i].strip():
                if res_numbers_int[j] == aa_position[i]:
                    scores.append(predictions2[i])
                    flag = False
                    # print(consurf_res_num[i], consurf_score_normalized[i])
                    # print(consurf_score_normalized[i])
                    print(predictions2[i])
            elif aa_format_dic[vars[j].rstrip()[:3]] != aa_seq[i].strip():
                if res_numbers_int[j] == aa_position[i]:
                    flag = False
                    print(' ')
        if flag == True:
            print(vars[j])
