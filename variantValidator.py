#! /usr/bin/env python3
import pandas as pd
import re
'''for variantValidator format input @ (https://batch.variantvalidator.org/batchvalidator/) to get the genomic coordinates for dbNSFP (Predictions) @ https://sites.google.com/site/jpopgen/dbNSFP'''

gene = 'cryba1', 'cryba2', 'cryba4', 'crybb1', 'crybb2', 'crybb3', 'crygc', 'crygd', 'crygs'
#######################################################
# f_vars = pd.read_excel('/Users/mdefsss2/other_genes/' + gene + '/' + gene + '_analysis.xlsx', sheet_name='analysis')
f_vars = pd.read_excel('/Users/mdefsss2/other_genes/crystallins/crystallins_analysis.xlsx')
output_file = open('/Users/mdefsss2/other_genes/crystallins/_variant_validator_input.txt', "w+")
file_contains_gene_NM = pd.read_excel('/Users/mdefsss2/other_genes/genes_list.xlsx')
hgvs = f_vars['HGVS']
variants = f_vars['variants']
genes = f_vars['genes']
three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}

'''extracting res number'''
# variant_residue_numbers = (re.findall(r'\d+', str(variants)))
# print('residue numbers', variant_residue_numbers)
# print('number of residues: ' + str(len(variant_residue_numbers)))

'''extracting mutated & mutants'''
mutated = []
mutants = []
for mut in variants:
    mutants.append(mut.rstrip()[-3:])
    mutated.append(mut.rstrip()[:3])
one_letter_mutated = []
one_letter_mutants = []
for each in range(len(mutants)):
    one_letter_mutated.append(three_letters_to_one_letter[mutated[each]])
    one_letter_mutants.append(three_letters_to_one_letter[mutants[each]])

for g in gene:
    refSeq = ''
    for i in range(len(file_contains_gene_NM)): #obtaining NM_number for the trascript used
        if file_contains_gene_NM['gene'][i] == g:
            refSeq += file_contains_gene_NM['refseq'][i]

    '''for transvar input (genomic coordinates needed to identify variants in dbSNFP)'''
    for i in range(len(variants)):
        if genes[i] == g:
            output_file.write(str(refSeq) + ':' + str(hgvs[i])+'\n')
            print(str(refSeq) + ':' + str(hgvs[i]))
output_file.close()
#######################################################
###########################################################
'''extra code'''
# f = open("/Users/mdefsss2/dbNSFP4/dbNSFP4.0a_variant.chrX", "r")
# f_coordinates = pd.read_csv('/Users/mdefsss2/revel/transvar_results.csv')
# coordinates = f_coordinates['coordinates(gDNA/cDNA/protein)']
# coords= []
# aa_pos= []
# for i in coordinates:
#     split_i = i.strip().split('/')
#     k= split_i[0]
#     l= split_i[1]
#     m= split_i[2]
#     coords.append(re.findall(r'\d+', str(k)))
#     aa_pos.append(re.findall(r'\d+', str(m)))
# coords_list = [''.join(x) for x in coords]
# aa_pos_list = [''.join(x) for x in aa_pos]
# # print('coordinates in a list: ', coords_list)
# # residue_numbers = data.columns.get_loc(variant_residue_numbers_list)
# coords_int = [int(i) for i in coords_list]
# aa_pos_int = [int(i) for i in aa_pos_list]
# # print('genomic coordinates: ', coords_int)
#
# # transcript = 'ENST00000376265'
# #####################################
# # output = open("/Users/mdefsss2/revel/revel_scores.txt", "w+")
# # counter = 0
# # # for j in range(len(coords_int)):
# # for line in f:
# #     cells = line.strip().split("\t")
# #     output.write(str(cells[8]) + ' ' + str(cells[12]) + ' ' + str(cells[4]) + ' ' + str(cells[5]) + "\n")
# #     print(str(cells[8]) + ' ' + str(cells[12]) + ' ' + str(cells[4]) + ' ' + str(cells[5]) + "\n")
# # output.close()
# # exit(True)
# ####################################
# output = open("/Users/mdefsss2/revel/cacna1f_subset.txt", "w+")
# counter = 0
#
# with open('/Users/mdefsss2/revel/revel_scores.txt', mode='r') as f:
#     data = f.read()
#     cells = data.splitlines()
#     for i in cells:
#         each = i.split()
#         for j in range(len(coords_int)):
#     # for line in f:
#         # if line.startswith("#chr"):
#         # cells = line.strip().split("\t")
#         # for i in cells:
#         # if str(cells[12]) == gene and str(cells[8]) == str(coords_int[j]):f
#             if str(each[0]) == str(coords_int[j]) and str(cells[1]) == gene and str(cells[2]) == one_letter_mutated[j] and str(cells[3]) == one_letter_mutants[j]:
#             # if str(cells[12]) == gene and str(cells[4]) == one_letter_mutated[j] and str(cells[5]) == one_letter_mutants[j]:
#
#             # if cells[12] == gene and cells[18] == str(hgvs[j]):
#             # if cells[8] or cells[10] == str(coords_int[j]):
#                 output.write(str(cells[0]) + ' ' + str(cells[1]) + ' ' + str(cells[2]) + ' ' + str(cells[3]) + "\n")
#                 counter+=1
#                 # print(str(cells[8]) + ' ' + str(cells[12]) + ' ' + str(cells[4]) + ' ' + str(cells[5]) + "\n")
#     output.close()
#     # break
# print('number of variants matched: ', counter)
# "-".join(cells[1], cells[3])