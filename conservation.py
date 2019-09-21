#!/usr/bin/env python

import pandas as pd
import re
from Bio import AlignIO
from PythonForBiology.MultipleAlign import profile
from PythonForBiology.Alignments import BLOSUM62
import numpy as np
'''21/04/18
creating a conservation score for residues in an alignment 
source used for the code: "cambridge core python programming for biology-chapters 12-14 : pairwise - multiple alignment"

the orthologues, returned by ENSEMBL for each gene, were used to calculate the conservation scores.
The orthologues from ALL SPECIES is downloaded (24/04/18) as fasta & unaligned with the query gene sequence added to the begining of the list. They were later aligned in SEAVIEW using clustal
However, since the code was not reading the clustal files and the alignment did not appear normal. THEREFORE, muscle was used (14/06) and conservation scores recorded in beta_gamma_hgmd_training_set.xlsx 

Ultimately, the orthologs provided by alamut (Deforche A., Blavier A. (2010). Systematic Building of Multiple Protein Alignments for Variant Interpretation Human Genome Meeting poster) was used for conservation. these were extracted and used here to measure conservation scores'''

'''CHANGE THE name of the gene & the dataset for different runs '''

# this_gene = 'single_gene_analysis'
# dataset = 'hgmd'
alignment_0 = AlignIO.read(open('/Users/mdefsss2/RDH5_variants/alamut_orthologs.aln'), "clustal") # to read the alignment
data = pd.read_excel('/Users/mdefsss2/RDH5_variants/rdh5_analysis_training_set.xlsx', index_col = [0], sheetname=1)
# written_file = '/Users/mdefsss2/RDH5_variants/rdh5_analysis_training_set.xlsx'
# genes = file_sheet1.columns.get_loc('genes')
data['conservation_alamut']= np.zeros(len(data))
index_conservation = data.columns.get_loc("conservation_alamut")
variants = data['variants'] # the column containing variants

residue_numbers = data.columns.get_loc('res_number')

# alignment_0 = AlignIO.read(open("/Users/mdefsss2/crystallins/cryba2_alamut_orth.fa"), "fasta") # to read one alignment

print("Alignment length %i" % alignment_0.get_alignment_length()) # alignment length
# print(alignment)
# print(alignment_0.format('clustal')) # to specify output format
alignment = []
for record in alignment_0 :
    alignment.append(record.seq)
#     print(record.seq + " " + record.id) # shows detailed info of all sequences (aligned)
print('the 1st sequence in alignment: \n ', alignment[0])

'''functions imported & parameters edited'''
def getConservation(align, simMatrix):
    conservation = []
    prof = profile(align)
    for compDict in prof:
      items = list(compDict.items())
      items.sort(key=lambda x: x[1])
      score = 0.0
      for resA, compA in items:
          for resB, compB in items:
              score += compA * compB * simMatrix[resA][resB]
      bestLetter = items[-1][0]
      maxScore = simMatrix[bestLetter][bestLetter]

      if maxScore == 0:
          score = 0
          # print('Score: ', score)
          conservation.append(score)
      else:
          score /= maxScore

          conservation.append(score)
    # print('Conservation length %i' % len(conservation))
    return conservation
getConservation(alignment, BLOSUM62)

def makeSimilarityString(align, simMatrix, thresholds):
    simString = []
    conservation = getConservation(align, simMatrix)
    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 = thresholds
    for score in conservation:
        if score >= t1:
            symbol = '1'
            # print(score)
        elif score >= t2 and score < t1:
            symbol = '2'
        elif score >= t3 and score < t2:
            symbol = '3'
        elif score >= t4 and score < t3:
            symbol = '4'
        elif score >= t5 and score < t4:
            symbol = '5'
        elif score >= t6 and score < t5:
            symbol = '6'
        elif score >= t7 and score < t6:
            symbol = '7'
        elif score >= t8 and score < t7:
            symbol = '8'
        # elif score >= t10 and score < t8:
        #     symbol = '9'
        else:
            symbol = '9'
        simString += symbol
    print('number of assigned scores: ' + str(len(simString)))
    return simString
symbols = makeSimilarityString(alignment, BLOSUM62, (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)) ##  calculating conservation scores
print('residue conservation scores assigned: ', len(symbols))

'''To extract the human seq from the multi seq alignment'''
for seq in alignment:
    # print(seq)
    # print(symbols)
    all_aa = seq[0:]
    break
print('number of residues in the sequence: ', len(all_aa))
seq_str = str(all_aa)

seq_res = []
for aa in seq_str:
    seq_res.append(aa) #  list of human residues from the alignment
print('the aligned residues & their numbers: ', seq_res, len(seq_res))

human_seq_res = []
for aa in seq_res:
    if aa != '-':
        human_seq_res.append(aa)
symbols_int = [int(float(i)) for i in symbols]
# print('the symbols: ', symbols_int)

'''combining the two lists of seq residues and their corresponding conservation scores'''
resNumber_residue_conScore_dict = {k: v for k, v in zip(enumerate(seq_res, start=1), symbols_int) for r in range(len(seq_res))}
# print('length of residue numbers & their conservation scores: ' + str(len(resNumber_residue_conScore_dict)))
# print('residue numbers & their conservation scores: ' + str(resNumber_residue_conScore_dict))

'''combining the two lists of seq residues and their corresponding conservation scores and removing the unmodelled residues through '-' '''
resNumber_residue_conScore_dict_human_seq = {k: v for k, v in zip(enumerate(seq_res, start=1), symbols_int) for r in range(len(seq_res)) if k[1] != '-'}
print('number of residues in the sequence & conservation scores assigned: ' + str(len(resNumber_residue_conScore_dict_human_seq)))
# print('HUMAN residue numbers & their conservation scores: ', resNumber_residue_conScore_dict_human_seq)

human_res_and_numbering_in_the_alignmt = list(key for key in resNumber_residue_conScore_dict if key[1] != '-') ##  only keeping the modelled residues but no conservation scores!
# print('length of modelled residues: ' + str(len(human_res_and_numbering_in_the_alignmt)))
# print('aligned modelled residues and cons scores: ' + str(human_res_and_numbering_in_the_alignmt))

'''getting the variants numbering through linking the alignment-sequence numbering with the variants '''

'''extracting res number'''
# variants = file_sheet1['variants']
# variant_residue_numbers = []
# for i in range(len(variants)):
#     gene = file_sheet1.iloc[i, genes]
#     if gene == this_gene:
#         res_number_for_each_var = file_sheet1.iloc[i, number_of_residues]
#         # variants_for_this_gene.append(res_number_for_each_var)
#         # print('variants: ', res_number_for_each_var)
#         variant_residue_numbers.append((re.findall(r'\d+', str(variants[i]))))
# print('the variants\' residue numbers', variant_residue_numbers)
# print('number of residues: ', len(variant_residue_numbers))

'''converting the above list of lists to a list of strings, then integers here'''
# variant_residue_numbers_list = [int(''.join(x)) for x in variant_residue_numbers]
# print('residue numbers: ', variant_residue_numbers_list)

'''constructing a list of the values which are the conservation scores'''
con_score = [v for k, v in resNumber_residue_conScore_dict_human_seq.items()]
print('conservation scores of the residues in the human sequence: ', con_score)

'''combining the two lists of residue numbers enumerated  and their corresponding conservation scores i.e. key is made up of a tuple'''
human_seq_resNumber_conScore_dict = {k: v for k, v in zip(enumerate(human_seq_res, start=1), con_score) for i in range(len(human_seq_res))}
# print('enumerated human residue numbers & conservation scores: ', human_seq_resNumber_conScore_dict)

'''constructing a list of the keys, which are made up of the human residue numbering and the residues'''
seq_res_num_con_score = [k for k, v in human_seq_resNumber_conScore_dict.items()]
print('the residues & their numbering in the sequence: ', seq_res_num_con_score)

'''extracting the residue numbers only'''
human_seq_res_num = [num[0] for num in seq_res_num_con_score]
# print('sequence res numbers: ' + str(human_seq_res_num))

'''combining the two lists of the human residue numbering and their corresponding conservation scores'''
human_residueNumber_conScore_dict = {k: v for k, v in zip(human_seq_res_num, con_score) for i in range(len(human_seq_res_num))}
print('human sequence residue numbers & their conservation scores: ', human_residueNumber_conScore_dict)

'''separating the residues of the sequence from the dictionary built earlier'''
human_seq_residues = [v for k, v in human_seq_resNumber_conScore_dict.keys()]
print('the residues in the human sequence: ', human_seq_residues)
print('number of residues: ', len(human_seq_residues))

variants_for_this_gene = []
for res in range(len(data)):
    # not_on_strand_flag= False
    # gene = data.iloc[res, genes]
    # if gene == this_gene:
    res_number_for_each_var = data.iloc[res, residue_numbers]
    variants_for_this_gene.append(res_number_for_each_var)
    # print('variants: ', res_number_for_each_var)

'''combining the variant list with the corresponding conservation scores '''
print('residue numbers: ', variants_for_this_gene)

print('the scores assigned to the variants\' residues: ')
variant_con_scores = []
for number in range(len(variants_for_this_gene)):
    score = human_residueNumber_conScore_dict[variants_for_this_gene[number]]
    variant_con_scores.append(score)
    # print(variants_for_this_gene[number])
    print(score)

'''writing the scores generated to the file [NOT CORRECT!!]'''
# for i in range(len(file_sheet1)):
#     gene = file_sheet1.iloc[i, genes]
#     if gene == this_gene:
#         for k in range(len(variant_con_scores)):
#             file_sheet1.iloc[k, index_conservation] = variant_con_scores[k]
# print('variant residues & their conservation scores: ', variants_residues_and_con_scores)

'''to write new values to the corresponding cell in a new column'''
# file_sheet1.to_excel(written_file, index=True) # re-writing the data

'''show all residues' conservation scores '''
import matplotlib.pyplot as plt
plt.bar(range(len(human_residueNumber_conScore_dict)),list(human_residueNumber_conScore_dict.values()), align='center')
plt.xticks(range(len(human_residueNumber_conScore_dict)), list(human_residueNumber_conScore_dict.keys()))
plt.show()
'''as a histogram'''
range=[1,2,3,4,5,6,7,8,9] #bin edges
plt.hist(con_score, rwidth=0.8, bins=range)
plt.show()