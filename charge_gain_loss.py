#! /usr/bin/env python3

import numpy as np
import pandas as pd
import re
import matplotlib as plt
import random

'''30/10/18 looking to find variants with change in charge, CHANGE 1) DATA'''

data = pd.read_excel('/Users/mdefsss2/CACNA1F/variants data/cacna1f_pathogenicMissense_vars.xlsx') # parameter 'index=False' skip automatic indices
variants_numbers = len(data)
print('length of the data: ', variants_numbers)
variants = data['Variants']
pd.options.display.max_rows=1000
three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}
R = 1
K = 2
H = 3
D = 4
E = 5
A = 6
N = 7
S = 8
Q = 9
T = 10
M = 11
G = 12
L = 13
I = 14
V = 15
F = 16
Y = 17
W = 18
C = 19
P = 20
# res_num_dic = {'Y': 1, 'I': 2, 'N': 3, 'T': 4, 'L': 5, 'Q': 6, 'C': 7, 'F': 8, 'W': 9, 'P':10, 'S': 11, 'A': 12, 'V': 13, 'G': 14, 'M': 15, 'D': 16, 'E': 17, 'K': 18, 'R': 19, 'H': 20}
positive_residues_list = ['K', 'R', 'H']
negative_residues_list = ['D', 'E']
uncharged_residues_list = ['Y', 'I', 'N', 'T', 'L', 'Q', 'C', 'F', 'W', 'P', 'S', 'A', 'V', 'G', 'M']
#################################################################### file 1:

# target_pathogenicity = data['Target_Patho']
# pathogenic_variants = []
# nonPathogenic_variants = []
# for i in range(len(variants)):
#     if target_pathogenicity[i] =='Pathogenic':
#         pathogenic_variants.append(variants[i])
#     if target_pathogenicity[i] =='NonPathogenic':
#         nonPathogenic_variants.append(variants[i])

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
observed_positive_gain_counter = 0
observed_negative_gain_counter = 0
observed_positive_loss_counter = 0
observed_negative_loss_counter = 0
observed_neutral_counter = 0
observed_other_counter = 0

for i in range(len(variants)):
    if one_letter_mutants_list[i] in positive_residues_list and one_letter_mutated_list[i] not in positive_residues_list:
        observed_positive_gain_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutants_list[i] in negative_residues_list and one_letter_mutated_list[i] not in negative_residues_list:
        observed_negative_gain_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutated_list[i] in positive_residues_list and one_letter_mutants_list[i] not in positive_residues_list:
        observed_positive_loss_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutated_list[i] in negative_residues_list and one_letter_mutants_list[i] not in negative_residues_list:
        observed_negative_loss_counter +=1
        # print(hgmd_training_data.index[i])
    elif one_letter_mutated_list[i] and one_letter_mutants_list[
        i] in uncharged_residues_list:
        observed_neutral_counter +=1
    else:
        observed_other_counter+=1
        # print(hgmd_training_data.index[i])
print('observed_positive_gain_counter: ', observed_positive_gain_counter)
print('observed_negative_gain_counter: ', observed_negative_gain_counter)
print('observed_positive_loss_counter: ', observed_positive_loss_counter)
print('observed_negative_loss_counter: ', observed_negative_loss_counter)
print('observed_neutral_counter: ', observed_neutral_counter)
print('observed_other_counter: ', observed_other_counter)

observed_positive_gain_fraction = float(observed_positive_gain_counter) / float(variants_numbers)
observed_negative_gain_fraction = float(observed_negative_gain_counter) / float(variants_numbers)
observed_positive_loss_fraction = float(observed_positive_loss_counter) / float(variants_numbers)
observed_negative_loss_fraction = float(observed_negative_loss_counter) / float(variants_numbers)

'''In regard to possible mutation outcomes:'''
total_residues = 20
neutral_residues = 15
negative_residues = 2
positive_residues = 3

'''Expected number of positively/negatively charged residues gained/lost & compare to those observed'''
expected_positive_residue_gain_fraction = observed_positive_gain_fraction * (neutral_residues + negative_residues) / total_residues
print ('Expected positive-gain: ', expected_positive_residue_gain_fraction)
print('observed positive-gain: ', observed_positive_gain_fraction)
expected_negative_residue_gain_fraction = observed_negative_gain_fraction * (neutral_residues + positive_residues) / total_residues
print ('Expected negative-gain: ', expected_negative_residue_gain_fraction)
print('observed negative-gain: ', observed_negative_gain_fraction)
expected_positive_residue_loss_fraction = observed_positive_loss_fraction * (neutral_residues + negative_residues) / total_residues
print ('Expected positive-loss: ', expected_positive_residue_loss_fraction)
print('observed positive-loss: ', observed_positive_loss_fraction)
expected_negative_residue_loss_fraction = observed_negative_loss_fraction * (neutral_residues + positive_residues) / total_residues
print ('Expected negative-loss: ', expected_negative_residue_loss_fraction)
print('observed negative-loss: ', observed_negative_loss_fraction)

'''Randomly generating residues & comparing it to observed mutant residues
 -counter_higher: the number of times that the generated fraction of gained charge mutations is higher than observed_fraction
 -counter_lower: the number of times that the generated fraction of gained charge mutations is less than or equal to observed_fraction'''
counter_positive_generated_higher_than_observed_gain = 0
counter_positive_generated_lower_than_observed_gain = 0
counter_negative_generated_higher_than_observed_gain = 0
counter_negative_generated_lower_than_observed_gain = 0
counter_positive_generated_higher_than_observed_loss = 0
counter_positive_generated_lower_than_observed_loss = 0
counter_negative_generated_higher_than_observed_loss = 0
counter_negative_generated_lower_than_observed_loss = 0
n_iterations = 10000
for j in range(n_iterations): # Repetitions allowed
    list_vals = []
    for m in range(len(variants)): # Generate residues/random numbers
        x = random.randint(1, 20)
        list_vals.append(x)
    neutral_counter = 0
    positive_counter = 0
    negative_counter = 0

    for val in list_vals:
        if val < 4:
            positive_counter += 1
        elif val ==4 or val ==5:
            negative_counter += 1
        else:
            neutral_counter += 1
    generated_positive_fraction = float(positive_counter) / float(positive_counter + negative_counter + neutral_counter)
    generated_negative_fraction = float(negative_counter) / float(positive_counter + negative_counter + neutral_counter)

    if generated_positive_fraction >= observed_positive_gain_fraction:
        counter_positive_generated_higher_than_observed_gain += 1
    elif generated_positive_fraction < observed_positive_gain_fraction:
        counter_positive_generated_lower_than_observed_gain += 1
    else:
       print("problem1")

    if generated_positive_fraction >= observed_positive_loss_fraction:
        counter_positive_generated_higher_than_observed_loss += 1
    elif generated_positive_fraction < observed_positive_loss_fraction:
        counter_positive_generated_lower_than_observed_loss += 1
    else:
       print("problem1")

    if generated_negative_fraction >= observed_negative_gain_fraction:
        counter_negative_generated_higher_than_observed_gain += 1
    elif generated_negative_fraction < observed_negative_gain_fraction:
        counter_negative_generated_lower_than_observed_gain += 1
    else:
       print("problem2")

    if generated_negative_fraction >= observed_negative_loss_fraction:
        counter_negative_generated_higher_than_observed_loss += 1
    elif generated_negative_fraction < observed_negative_loss_fraction:
        counter_negative_generated_lower_than_observed_loss += 1
    else:
       print("problem2")

p_value_positive_gain = float(counter_positive_generated_higher_than_observed_gain) / float(counter_positive_generated_lower_than_observed_gain + counter_positive_generated_higher_than_observed_gain)
p_value_negative_gain = float(counter_negative_generated_higher_than_observed_gain) / float(counter_negative_generated_lower_than_observed_gain + counter_negative_generated_higher_than_observed_gain)
print ("Positive gain p-value:", p_value_positive_gain)
print ("Negative gain p-value:", p_value_negative_gain)

p_value_positive_loss = float(counter_positive_generated_higher_than_observed_loss) / float(counter_positive_generated_lower_than_observed_loss + counter_positive_generated_higher_than_observed_loss)
p_value_negative_loss = float(counter_negative_generated_higher_than_observed_loss) / float(counter_negative_generated_lower_than_observed_loss + counter_negative_generated_higher_than_observed_loss)
print ("Positive loss p-value:", p_value_positive_loss)
print ("Negative loss p-value:", p_value_negative_loss)

####################################################################
