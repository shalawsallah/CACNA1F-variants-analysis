#! /usr/bin/env python3
import pandas as pd
import random
'''14/04/18 change in charge amongst core or surface residues'''

data = pd.read_excel('/Users/mdefsss2/cacna1f/cacna1f_analysis.xlsx') # parameter 'index=False' skips automatic indices
all_analysis = False
surface_only_analysis = True
core_only_analysis = False
surface_core_threshold = 10
###############################################
solvent_access = data['rel_side_acc'].tolist()
naccess_scores = data.columns.get_loc('rel_side_acc')
variants = data['variants'].tolist()

positive_residues_list = ['K', 'R', 'H']
negative_residues_list = ['D', 'E']
uncharged_residues_list = ['Y', 'I', 'N', 'T', 'L', 'Q', 'C', 'F', 'W', 'P', 'S', 'A', 'V', 'G', 'M']

'''In regard to possible mutation outcomes'''
total_residues = 20
neutral_residues = 15
negative_residues = 2
positive_residues = 3

three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
                               'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
                               'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}

'''extracting mutants & mutated for variants with a low solvent accessibility score'''
mutant_residues_list = []
mutated_residues_list = []

for mut in range(len(variants)):
    h_naccess = data.iloc[mut, naccess_scores]
    if surface_only_analysis == True:
        if h_naccess >= surface_core_threshold:
            mutant_residues_list.append(variants[mut].rstrip()[-3:])
            mutated_residues_list.append(variants[mut].rstrip()[:3])
    elif core_only_analysis == True:
        if h_naccess < surface_core_threshold:
            mutant_residues_list.append(variants[mut].rstrip()[-3:])
            mutated_residues_list.append(variants[mut].rstrip()[:3])
    elif all_analysis == True:
        mutant_residues_list.append(variants[mut].rstrip()[-3:])
        mutated_residues_list.append(variants[mut].rstrip()[:3])
print('number of mutated: ', len(mutated_residues_list))
print('number of mutants: ', len(mutant_residues_list))
# print('mutants: ', mutant_residues_list)

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
observed_charge_flip_counter = 0
observed_neutral_counter = 0
# other_counter = 0

for i in range(len(mutated_residues_list)):
    if one_letter_mutants_list[i] in positive_residues_list and one_letter_mutated_list[i] in uncharged_residues_list:
        observed_positive_gain_counter +=1
        print('positive gain')
    elif one_letter_mutated_list[i] in positive_residues_list and one_letter_mutants_list[i] in uncharged_residues_list:
        observed_positive_loss_counter +=1
        print('positive loss')
    elif one_letter_mutants_list[i] in negative_residues_list and one_letter_mutated_list[i] in uncharged_residues_list:
        observed_negative_gain_counter +=1
        print('negative gain')
    elif one_letter_mutated_list[i] in negative_residues_list and one_letter_mutants_list[i] in uncharged_residues_list:
        observed_negative_loss_counter +=1
        print('negative loss')
    elif (one_letter_mutated_list[i] in negative_residues_list and one_letter_mutants_list[i] in positive_residues_list) or (one_letter_mutated_list[i] in positive_residues_list and one_letter_mutants_list[i] in negative_residues_list):
        observed_charge_flip_counter +=1
        print('charge flip')
    elif one_letter_mutated_list[i] and one_letter_mutants_list[
        i] in uncharged_residues_list or one_letter_mutated_list[i] and one_letter_mutants_list[
        i] in negative_residues_list or one_letter_mutated_list[i] and one_letter_mutants_list[
        i] in positive_residues_list:
        observed_neutral_counter +=1
        print('none')
    else:
        print('?? mutation missed out!')

print('######################################')
print('observed_positive_gain_counter: ', observed_positive_gain_counter)
print('observed_negative_gain_counter: ', observed_negative_gain_counter)
print('observed_positive_loss_counter: ', observed_positive_loss_counter)
print('observed_negative_loss_counter: ', observed_negative_loss_counter)
print('observed_charge_flip_counter: ', observed_charge_flip_counter)
print('observed_neutral_counter: ', observed_neutral_counter)
all_counters = observed_negative_gain_counter+observed_negative_loss_counter+observed_neutral_counter+observed_positive_gain_counter+observed_positive_loss_counter+observed_charge_flip_counter
print('combined counters: ', all_counters)
print('number of variants: ', len(mutated_residues_list))

observed_positive_gain_fraction = float(observed_positive_gain_counter) / float(len(mutated_residues_list))
observed_negative_gain_fraction = float(observed_negative_gain_counter) / float(len(mutated_residues_list))
observed_positive_loss_fraction = float(observed_positive_loss_counter) / float(len(mutated_residues_list))
observed_negative_loss_fraction = float(observed_negative_loss_counter) / float(len(mutated_residues_list))
observed_charge_flip_fraction = float(observed_charge_flip_counter) / float(len(mutated_residues_list))

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
