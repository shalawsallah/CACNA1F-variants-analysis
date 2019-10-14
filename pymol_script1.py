#! /usr/bin/env python3

import pandas as pd
import re
import xlrd

'''29/10/18 residue size change calculation for variants

PYMOL SCRIPT 1. '''

file = pd.read_csv('/Users/mdefsss2/Desktop/weka_data/combined_train_test_data.csv', index_col = [0])
pd.options.display.max_rows=1000 # to display, and later use as input, all of the rows, and not just a subset that pandas automatically generates

variants = file.index
class_target = file['Target_Patho']

'''extracting res number'''
# variant_residue_numbers = (re.findall(r'\d+', str(variants)))
# print('residue numbers', variant_residue_numbers)
# print('number of residues: ' + str(len(variant_residue_numbers)))

pathogenic_variant_residue_numbers = []
nonpathogenic_variant_residue_numbers = []
for i in range(len(variants)):
    if class_target[i] == 'Pathogenic':
        res_numbers = (re.findall(r'\d+', str(variants[i])))
        pathogenic_variant_residue_numbers.append(res_numbers)
    elif class_target[i] == 'NonPathogenic':
        res_numbers = (re.findall(r'\d+', str(variants[i])))
        nonpathogenic_variant_residue_numbers.append(res_numbers)
pathogenic_variant_residue_numbers_list = [''.join(x) for x in pathogenic_variant_residue_numbers]
nonpathogenic_variant_residue_numbers_list = [''.join(x) for x in nonpathogenic_variant_residue_numbers]

for res in pathogenic_variant_residue_numbers_list:
    print(int(res))
print('number of pathogenic residues: ', len(pathogenic_variant_residue_numbers_list))

for res in nonpathogenic_variant_residue_numbers_list:
    print(int(res))
print('number of non-pathogenic residues: ', len(nonpathogenic_variant_residue_numbers_list))

'''the results from here are stored in a txt file for Script 2'''
