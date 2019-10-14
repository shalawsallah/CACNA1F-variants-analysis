#! /usr/bin/env python3

import numpy as np
import pandas as pd
import re
import matplotlib as plt
import random

'''30/10/18 looking to find variants resulting in hydrophobicity changes.
CHANGE the file only, for new data'''

data = pd.read_excel('/Users/mdefsss2/CACNA1F/variants data/cacna1f_pathogenicMissense_vars.xlsx') # parameter 'index=False' skip automatic indices
variants_numbers = len(data)
print('length of the data: ', variants_numbers)
variants = data['Variants']
pd.options.display.max_rows=1000
variants_numbers = len(data)
print('length of the data: ', variants_numbers)

three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}
# res_num_dic = {'Y': 1, 'I': 2, 'N': 3, 'T': 4, 'L': 5, 'Q': 6, 'C': 7, 'F': 8, 'W': 9, 'P':10, 'S': 11, 'A': 12, 'V': 13, 'G': 14, 'M': 15, 'D': 16, 'E': 17, 'K': 18, 'R': 19, 'H': 20}
non_hydrophobic_residues_list = ['K', 'R', 'H', 'D', 'E', 'Y', 'N', 'Q', 'P', 'T', 'S', 'A', 'G']
hydrophobic_residues_list = ['I', 'L', 'V', 'M', 'F', 'W', 'C'] # replacing (Steitz and Goldman, CACNA1F with M.J. Betts, R.B. Russell. Amino acid properties and consequences of subsitutions. In Bioinformatics for Geneticists, M.R. Barnes, I.C. Gray eds, Wiley, 2003

'''In regard to possible mutation outcomes'''
total_residues = 20
non_hydrophobic_residues = 13
hydrophobic_residues = 7
####################################################################
# residues numbered 1-7 are hydrophobic :

M = 1
L = 2
I = 3
V = 4
F = 5
W = 6
C = 7
A = 8
S = 9
T = 10
G = 11
Y = 12
P = 13
Q = 14
N = 15
R = 16
K = 17
H = 18
D = 19
E = 20

'''extracting res number'''
# variant_residue_numbers = (re.findall(r'\d+', str(data.index)))
# print('residue numbers', variant_residue_numbers)
# print('number of residues: ' + str(len(variant_residue_numbers)))
'''extracting mutants'''
mutant_residues_list = []
for mut in variants:
    mutant_residues_list.append(mut.rstrip()[-3:])
print('number of mutants: ', len(mutant_residues_list))
# print('mutants: ', mutant_residues_list)
'''extracting mutated'''
mutated_residues_list = []
for mutd in variants:
    mutated_residues_list.append(mutd[:3])
print('number of mutated: ', len(mutated_residues_list))
one_letter_mutants_list = []
for each in range(len(mutant_residues_list)):
    one_letter_mutants_list.append(three_letters_to_one_letter[mutant_residues_list[each]])
# print('mutants: ', one_letter_mutants_list)
one_letter_mutated_list = []
for each in range(len(mutated_residues_list)):
    one_letter_mutated_list.append(three_letters_to_one_letter[mutated_residues_list[each]])
# print('mutated: ', one_letter_mutated_list)

'''calculating the number of positively/negatively charged residues gained/lost from the list of variants'''

observed_hydrophob_gain_counter = 0
observed_hydrophob_loss_counter = 0
observed_neutral_counter = 0
for i in range(len(variants)):
    if one_letter_mutants_list[i] in hydrophobic_residues_list and one_letter_mutated_list[i] not in hydrophobic_residues_list:
        observed_hydrophob_gain_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutated_list[i] in hydrophobic_residues_list and one_letter_mutants_list[i] not in hydrophobic_residues_list:
        observed_hydrophob_loss_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutated_list[i] in hydrophobic_residues_list and one_letter_mutants_list[
        i] in hydrophobic_residues_list:
        observed_neutral_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutants_list[i] in non_hydrophobic_residues_list and one_letter_mutated_list[i] in non_hydrophobic_residues_list:
        observed_neutral_counter +=1
        # print(hgmd_training_data.index[i])
print('neutral counter: ', observed_neutral_counter, 'gain counter: ', observed_hydrophob_gain_counter, 'loss counter: ', observed_hydrophob_loss_counter)

observed_hydrophobic_gain_fraction = float(observed_hydrophob_gain_counter) / float(variants_numbers)
observed_hydrophobic_loss_fraction = float(observed_hydrophob_loss_counter) / float(variants_numbers)


'''Expected number of hydrophobic residues gained/lost compared to those observed'''
expected_hydrophobic_residue_gain_fraction = observed_hydrophobic_gain_fraction * (non_hydrophobic_residues) / total_residues
print('Expected hydrophobic gain fraction: ', expected_hydrophobic_residue_gain_fraction)
print("Observed hydrophobic gain fraction:", observed_hydrophobic_gain_fraction)

expected_hydrophobic_residue_loss_fraction = observed_hydrophobic_loss_fraction * (hydrophobic_residues) / total_residues
print('Expected hydrophobic loss fraction: ', expected_hydrophobic_residue_loss_fraction)
print("Observed hydrophobic loss fraction:", observed_hydrophobic_loss_fraction)

'''Randomly generating residues & comparing it to observed mutant residues
- counter_higher: the number of times that the generated fraction of gained hydrophobic residues is higher than observed_fraction
- counter_lower: the number of times that the generated fraction of gained hydrophobic residues is less than or equal to observed_fraction'''
counter_hydrophobic_generated_higher_than_observed_gain = 0
counter_hydrophobic_generated_lower_than_observed_gain = 0
counter_non_hydrophobic_generated_higher_than_observed_loss = 0
counter_non_hydrophobic_generated_lower_than_observed_loss = 0
n_iterations = 10000

for j in range(n_iterations): # Repetitions allowed
    list_vals = []
    for m in range(len(variants)):     # Generate random numbers/residues
        x = random.randint(1, 20)
        list_vals.append(x)
    non_hydrophobic_counter = 0
    hydrophobic_counter = 0

    for val in list_vals:
        if val > 7:
            non_hydrophobic_counter += 1
        else:
            hydrophobic_counter += 1

    generated_hydrophobic_fraction = float(hydrophobic_counter) / float(hydrophobic_counter + non_hydrophobic_counter)
    generated_non_hydrophobic_fraction = float(non_hydrophobic_counter) / float(hydrophobic_counter + non_hydrophobic_counter)

    if generated_hydrophobic_fraction >= observed_hydrophobic_gain_fraction:
        counter_hydrophobic_generated_higher_than_observed_gain += 1
    elif generated_hydrophobic_fraction < observed_hydrophobic_gain_fraction:
        counter_hydrophobic_generated_lower_than_observed_gain += 1
    else:
       print("problem1")

    if generated_non_hydrophobic_fraction >= observed_hydrophobic_loss_fraction:
        counter_non_hydrophobic_generated_higher_than_observed_loss += 1
    elif generated_non_hydrophobic_fraction < observed_hydrophobic_loss_fraction:
        counter_non_hydrophobic_generated_lower_than_observed_loss += 1
    else:
       print("problem2")

p_value_hydrophobic_gain = float(counter_hydrophobic_generated_higher_than_observed_gain) / float(counter_hydrophobic_generated_lower_than_observed_gain + counter_hydrophobic_generated_higher_than_observed_gain)
p_value_hydrophobic_loss = float(counter_non_hydrophobic_generated_higher_than_observed_loss) / float(counter_non_hydrophobic_generated_lower_than_observed_loss + counter_non_hydrophobic_generated_higher_than_observed_loss)

print("hydrophobic gain p-value:", p_value_hydrophobic_gain)
print("hydrophobic loss p-value:", p_value_hydrophobic_loss)
