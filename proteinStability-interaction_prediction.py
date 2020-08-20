#! /usr/bin/env python3

import pandas as pd
import numpy as np
import re
# import matplotlib as plt
import random

'''24/04/19 
extracting variants in the format needed to submit to mCSM @ http://biosig.unimelb.edu.au/mcsm/protein_protein
to obtain predictions on effect of mutations on protein folding/interaction energy

to ADD CHAIN ID to the PDB if one is missing: http://www.canoz.com/sdh/renamepdbchain.pl OR using CHIMERA'''
gene = 'rpgr'
template = '4jhn'
pdb = open('/Users/mdefsss2/x_linked_genes/'+gene+'/'+gene+'-model_'+template+'.pdb')
variant_file = pd.read_excel('/Users/mdefsss2/x_linked_genes/'+gene+'/'+gene+'_analysis.xlsx')
chain_id = 'A' #C for bb3
#####################################################
# for line in pdb: #info from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
#     chain_id.join(line[21].strip())
#     break
residues_in_struc = []
for line in pdb: #info from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    if line.startswith('ATOM'):
        res_num = int(line[22:26].strip())
        chainID = line[21:22].strip()
        if chainID == chain_id:
            residues_in_struc.append(res_num)
variants = variant_file['variants']
# genes= variant_file['genes']
modelled_res = variant_file['modelled']

'''extracting res number'''
variant_residue_numbers = []
for i in range(len(variants)):
    variant_residue_numbers.append((re.findall(r'\d+', str(variants[i]))))
# # print('residue numbers', variant_residue_numbers)
# print('number of residues: ', len(variant_residue_numbers))
#
# '''converting the above list of lists to a list of strings here'''
variant_residue_numbers_list = [''.join(x) for x in variant_residue_numbers]
# print('residue numbers: ', variant_residue_numbers_list)
#
# '''extracting mutants'''
# mutants = []
# for mut in variants:
#     mutants.append(mut.rstrip()[-3:])
# print(len(mutants))

three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}
for i in range(len(variants)):
    variant_residue_numbers = variant_residue_numbers_list[i]
    mutated = three_letters_to_one_letter[variants[i].rstrip()[:3]]
    mutants = three_letters_to_one_letter[variants[i].rstrip()[-3:]]
    modelled = modelled_res[i]
    if int(variant_residue_numbers) in residues_in_struc:
        print(chain_id, str(mutated) + str(variant_residue_numbers) + str(mutants))
    # elif int(variant_residue_numbers) not in residues_in_struc:
    #     print(' ')
# for line in pdb: #info from https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
#         chain_id = line[21].strip()
#         # print(chain_id)
#         if chain_id == 'A' or chain_id == 'X':
#             strand_start = int(line[23:26].strip())
#             # if renumbered == True:
#             #     strand_start = int(line[23:26].strip()) +1 #renummbered by +1
#             strand_start_res_name = line[17:20].strip()
#             strand_number_in_sheet = line[8:10].strip() #in each sheet
#             sheet_identifier = line[12:14].strip()
#             number_of_strands_in_sheet = line[15:16].strip()
#             strand_end = int(line[34:37].strip())
#             strand_end_res_name = line[28:31].strip()
#
#     if line[:5].strip() == 'HELIX':
#         # print(line)
#         chain_id = line[19].strip()
#         # print(chain_id)
#         if chain_id == 'A' or chain_id == 'X':
#             helix_start = int(line[22:25].strip())
#             # if renumbered == True:
#             #     helix_start = int(line[22:25].strip()) +1 #renummbered by +1
#             helix_start_res_name = line[15:18].strip()
#             helix_number = line[8:10].strip() #in each sheet
#             helix_end = int(line[34:37].strip())
#             # if renumbered == True:
#             #     helix_end = int(line[34:37].strip()) +1 #renummbered by +1
#             helix_end_res_name = line[27:30].strip()
#             helix_starts.append(helix_start)
#             helix_ends.append(helix_end)
#             # print(helix_start,helix_start_res_name, helix_end, helix_end_res_name)
# print('this gene: ', this_gene)
# print('start of strands: ', strand_starts)
# print('end of strands: ', strand_ends)
# print('start of helices: ', helix_starts)
# print('end of helices: ', helix_ends)
#
