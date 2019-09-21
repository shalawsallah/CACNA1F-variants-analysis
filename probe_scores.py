#! /usr/bin/env python3
import sys

import pandas as pd
import re
import os
from subprocess import call

'''06/04/18
Probe score calculation '''
##################################################################################
'''CHANGE /run_probe command, i.e. if ONE or MULTIPLE CHAINS are present & change the next 4 variables for each new file'''

data = pd.read_excel('/Users/mdefsss2/gpr143/gpr143_analysis.xlsx')
# gene= 'ndp'
# pdb_to_be_reduced = 'ocrl_model.pdb'
reduced_file_name = 'gpr143_raptorX_reduced.pdb'
new_gene_directory_name = 'mutrot_files_gpr143_vars'
################################################################
pd.options.display.max_rows=1000 # to display, and later use as input, all of the rows, and not just a subset that pandas automatically generates
variants = data['variants'].tolist() # the column containing variants
print(len(variants))
# print(variants)

'''extracting res number'''
variant_residue_numbers = []
for i in range(len(variants)):
    variant_residue_numbers.append((re.findall(r'\d+', str(variants[i]))))
# print('residue numbers', variant_residue_numbers)
print('number of residues: ', len(variant_residue_numbers))

'''converting the above list of lists to a list of strings here'''
variant_residue_numbers_list = [''.join(x) for x in variant_residue_numbers]
print('residue numbers: ', variant_residue_numbers_list)

'''extracting mutants'''
mutants = []
for mut in variants:
    mutants.append(mut.rstrip()[-3:])
print('number of mutants: ', mutants)
print(len(mutants))

three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}

'''Calling Mutrot and Probe and getting highest probe score for the variants'''
# os.system('cd /Users/mdefsss2/Desktop/test_mutrot_dir/') # to change directory???

'''the next 3 steps can also be done manually'''
'''to change directory into ./gene where pdb files are '''
# if __name__ == '__main__':
#     os.chdir("/Users/mdefsss2/" + gene)
#     os.system("pwd")
#     print ('dir changed')

'''adding H to the pdb file before analysis '''
# reduce_str = str('reduce -build ' + pdb_to_be_reduced + ' > ' + reduced_file_name)
# print(reduce_str)
# call(reduce_str, shell=True)

'''moving the reduced file to /test_mutrot directory'''
# move_str = str('mv ' + reduced_file_name + ' /Users/mdefsss2/Desktop/test_mutrot_dir/')
# print(move_str)
# call(move_str, shell=True)

'''to change directory into /test_mutrot where pdb & other files are '''
# if __name__ == '__main__':
os.chdir("/Users/mdefsss2/Desktop/test_mutrot_dir/")
os.system("pwd")
print ('dir changed')
# os.system("/bin/bash") #to remain in the directory after exiting the script or program

'''call and system are similar commands'''
os.system('mkdir ' + new_gene_directory_name)

for each in range(len(variant_residue_numbers_list)):
    # print('virtual memory: ', psutil.virtual_memory())
    # print('cpu percent: ', psutil.cpu_percent())
    mutrot_str = str('/Users/mdefsss2/Desktop/test_mutrot_dir/mutrot ' + reduced_file_name + ' -res ' + variant_residue_numbers_list[each] + ' -to ' + three_letters_to_one_letter[mutants[each]])
    # print(mutrot_str)
    call(mutrot_str, shell=True)
    probe_str = str('/Users/mdefsss2/Desktop/test_mutrot_dir/run_probe ' + variant_residue_numbers_list[each] + ' pdb_rot_out* | awk \'$0~\"grand tot\"||$0~\"^pdb_rot\" {print}\' -  > scores' + str(each) + '.txt')
    print(probe_str)
    call(probe_str, shell=True)
    call('mv pdb_rot_out* ' + new_gene_directory_name, shell=True)

for each in range(len(variant_residue_numbers_list)):
    highest_str = ('/Users/mdefsss2/Desktop/test_mutrot_dir/get_highest scores' + str(each) + '.txt')
    # print(highest_str)
    call(highest_str, shell=True)

call('mv scores* ' + new_gene_directory_name, shell=True)
call('mv rotamers.txt ' + new_gene_directory_name, shell=True)

