#! /usr/bin/env python3

import pandas as pd
import matplotlib as plt
from scipy import stats

'''06/04/18 residue size change calculation as a result of mutation'''

file = pd.read_excel('/Users/mdefsss2/cacna1f/cacna1f_analysis.xlsx')
variants = file['variants']
three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N', 'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R', 'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}
aa_volumes = {'G': 66, 'A': 92, 'S': 99, 'C': 106, 'T': 122, 'D': 125, 'P': 129, 'N': 135, 'V': 142, 'E': 155, 'Q': 161, 'H': 167, 'L': 168, 'I': 169, 'M': 171, 'K': 171, 'F': 203, 'Y': 203, 'R': 225, 'W': 240}

'''extracting mutated & mutants'''
mutated = []
mutants = []
for mut in variants:
    mutants.append(mut.rstrip()[-3:])
    mutated.append(mut.rstrip()[:3])
print('number of mutants: ', len(mutants))
print('number of mutated: ' , len(mutated))

one_letter_mutants = []
for each in range(len(mutants)):
    one_letter_mutants.append(three_letters_to_one_letter[mutants[each]])
print('one letter mutants: ', one_letter_mutants)

one_letter_mutated = []
for each in range(len(mutated)):
    one_letter_mutated.append(three_letters_to_one_letter[mutated[each]])
print('one letter mutated: ', one_letter_mutated)

for i in range(len(mutated)):
    sidechain_vol_difference = (aa_volumes[one_letter_mutated[i]] - aa_volumes[one_letter_mutants[i]])
    # print(mutated[i])
    print(sidechain_vol_difference)
#####################################################################################################################





