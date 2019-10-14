#! /usr/bin/env python3
import pandas as pd
import re
import xlrd
'''PYMOL SCRIPT 1. (& 2) '''
'''extracting the residue numbers from the variants to use in pymol

CHANGE:
-file
-pathogenic_or_benign class'''

file = pd.read_excel('/Users/mdefsss2/cacna1f/cacna1f_analysis.xlsx', sheet_name=0)
pathogenic_or_benign = 'pathogenic'
# pd.options.display.max_rows=1000 # to display, and later use as input, all of the rows, and not just a subset that pandas automatically generates
variants = file['variants']
pathogenicity = file['class']
pathogenicOrBenign_vars = []
for i in range(len(pathogenicity)):
    if pathogenicity[i] == pathogenic_or_benign:
        pathogenicOrBenign_vars.append(variants[i])

'''extracting res number'''
# variant_residue_numbers = (re.findall(r'\d+', str(variants)))
# print('residue numbers', variant_residue_numbers)
# print('number of residues: ' + str(len(variant_residue_numbers)))
variant_residue_numbers = []
for i in range(len(pathogenicOrBenign_vars)):
    res_numbers = (re.findall(r'\d+', str(pathogenicOrBenign_vars[i])))
    # print(res_numbers)
    variant_residue_numbers.append(res_numbers)
variant_residue_numbers_list = [''.join(x) for x in variant_residue_numbers]
# print(variant_residue_numbers_list)

res_list = []
print('number of residues: ', len(variant_residue_numbers_list))
for res in variant_residue_numbers_list:
    # print(int(res))
    res_list.append(int(res))
print('with duplicated residues mutated: ', len(res_list))

'''removing more that one mutation in the same residues (pymol crashes when giving it many residues) '''
dup_removed = []
for i in range(len(res_list)):
    if i > 0:
        if res_list[i] != res_list[i-1]:
            dup_removed.append(res_list[i])
            # print(res_list[i])
    # i = 0
    else:
        dup_removed.append(res_list[i])

print('following duplicated residues being removed: ', len(dup_removed))

array = ''
for i in range(len(dup_removed)):
    x= (str(dup_removed[i]) + '+')
    array = array + x
print(array[:-1])
