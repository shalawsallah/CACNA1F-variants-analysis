#!/usr/bin/env python

'''created 10/05/17, edited last 15/05/17
to ultimately measure a conservation score for the 24 orthologues of cacna1f.
code used here, is taken from "cambridge core python programming for biology-chapters 12-14 : pairwise - multiple alignment"'''

# from Bio.Alphabet import IUPAC
# iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"]
#
# def read_clustal_alignment(filename):
#     """ Read in the alignment stored in the CLUSTAL file, filename. Return
#     two lists: the names and sequences. """
#     names = []
#     alignment = []
#     f = open(filename)
#     for line in f:
#         line = line[:-1]
#         if len(line) == 0: continue
#         if '*' in line: continue
#         if 'CLUSTAL' in line: continue
#         t = line.split()
#         if len(t) == 2 and t[1][0] in iupac_alphabet:
#             if t[0] not in names:
#                 names.append(t[0])
#                 alignment.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
#             else:
#                 alignment[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X',
#                                                                                                          '-').replace(
#                     '\r', '')
#     return names, alignment
#
# alignment = read_clustal_alignment("/Users/mdefsss2/CACNA1F/24_orthologues/24_cacna1f_clustal_aligned.aln")
# print(alignment)

from Bio import AlignIO
from PythonForBiology.MultipleAlign import profile
alignment_0 = AlignIO.read(open("/Users/mdefsss2/CACNA1F/24_orthologues/24_cacna1f_clustal_aligned.aln"), "clustal") # to read one alignment
print("Alignment length %i" % alignment_0.get_alignment_length()) # alignment length
# print(alignment)
# print(alignment_0.format('clustal')) # to specify output format
alignment = []
for record in alignment_0 :
    alignment.append(record.seq)
    # print(record.seq + " " + record.id) # shows detailed info of all sequences (aligned)
# print('alignment: ' + str(alignment))

# we can use either consensus of profile for alignment
# def consensus(alignment, threshold=0.15):
#     n = len(alignment[0])
#     nSeq = float(len(alignment))
#     consensus = ''
#     for i in range(n):
#         counts = {}
#         for seq in alignment:
#             letter = seq[i]
#             if letter == '-':
#                 continue
#             counts[letter] = counts.get(letter, 0) + 1
#         fractions = []
#         for letter in counts:
#             frac = counts[letter] / nSeq
#             fractions.append([frac, letter])
#         fractions.sort()
#         bestFraction, bestLetter = fractions[-1]
#         if bestFraction < threshold:
#             consensus += 'X'
#         else:
#             consensus += bestLetter
#     return consensus

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
      score /= maxScore
      conservation.append(score)
    # print('Conservation length %i' % len(conservation))
    return conservation

from PythonForBiology.Alignments import BLOSUM62

getConservation(alignment, BLOSUM62)


def makeSimilarityString(align, simMatrix, thresholds):
    simString = []
    conservation = getConservation(align, simMatrix)
    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 = thresholds
    for score in conservation:
        # if score >= t1:
        #     symbol = '*'
        # elif score >= t2:
        #     symbol = ':'
        # elif score >= t3:
        #     symbol = '.'
        if score >= t1:
            symbol = '1'
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
    print('number of aa: ' + str(len(simString)))
    return simString
symbols = makeSimilarityString(alignment, BLOSUM62, (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)) ##  calculating conservation scores

# residue_number = 74 # where n = residue number
# mutated_residue = residue_number-1
# for seq in alignment:
#     print('start ' + seq + 'end')
#     print('start ' + symbols + 'end')
    # for res in seq:
    #     print(res)

    #print(type(seq))measureConservationScore.symbols
    #myseq=str(seq)#.split('\n')#replace('\n', '')
    #print(type(myseq))
    #for res in seq:
        #print(type(res))
        #if res != "\n":
            #myseq += res
    #print( myseq )

    #print (myseq)
    #print ("\n")
    #print (myseq.replace('A', 'X'))
    #break

# print(len(seq))
# print(len(symbols))
# variant_number = 164
############################################################################################################################

for seq in alignment:##  To extract the human cacna1f seq from the multi seq alignment

    # print(seq)
    # print(symbols)
    all_aa = seq[0:]
    break
# print(str(len(all_aa)))
seq_str = str(all_aa)

seq_res = []
for aa in seq_str:
    seq_res.append(aa) #  list of human residues from the alignment
# print(seq_res)
symbols_int = [float(i) for i in symbols]
# print(symbols_int)

# resNumber_residue_conScore_dict = {k: v for k, v in zip(enumerate(seq_res, start=1), symbols_int) for r in range(len(seq_res))} ##  combining the two lists of seq residues and their corresponding conservation scores

resNumber_residue_conScore_dict = {k: v for k, v in zip(enumerate(seq_res, start=1), symbols_int) for r in range(len(seq_res)) if k[1] != '-'} ##  combining the two lists of seq residues and their corresponding conservation scores and removing the unmodelled residues through '-'
print('length of residue numbers & their conservation scores: ' + str(len(resNumber_residue_conScore_dict)))
print('residue numbers & their conservation scores: ', '\n', resNumber_residue_conScore_dict)

# modelled_res = list(key for key in resNumber_residue_conScore_dict if key[1] != '-') ##  only keeping the modelled residues but no conservation scores!
# print('length of modelled residues: ' + str(len(modelled_res)))
#
# # print(type(modelled_res))
# print('aligned modelled residues and cons scores: ' + str(modelled_res))


############################################################################################################################

    # print(alignment[0])
# print('number of sequences: ' + str(len(alignment)))
# print('length of sequences: ' + str(len(seq)))
# print('residue: ' + seq[variant_number])
# print('cons score: ' + symbols[variant_number])
##############################################################

##  To differentiate between types of amino acids:

# AA_CATEGORIES = [('G','G'), ('A','A'), ('I','I'), ('V','V'),
#                  ('L','L'), ('M','M'), ('F','F'), ('Y','Y'),
#                  ('W','W'), ('H','H'), ('C','C'), ('P','P'),
#                  ('K','K'), ('R','R'), ('D','D'), ('E','E'),
#                  ('Q','Q'), ('N','N'), ('S','S'), ('T','T'),
#                  ('-','-'),
#                  ('acidic',    'DE'),
#                  ('hydroxyl',  'ST'),
#                  ('aliphatic', 'VIL'),
#                  ('basic', 'KHR'),
#                  ('tiny', 'GAS'),
#                  ('aromatic', 'FYWH'),
#                  ('charged', 'KHRDE'),
#                  ('small', 'AGSVTDNPC'),
#                  ('polar', 'KHRDEQNSTC'),
#                  ('hydrophobic','IVLFYWHAGMC'),
#                  ('turnlike',   'GASHKRDEQNSTC'),
#                  ('undef',      'AGVILMFYWHCPKRDEQNST-')]
#
# def getAlignProperties(align, categories):
#   properties = []
#   prof = profile(align)
#   for fracDict in prof:
#       letters = fracDict.keys()
#       for name, group in categories:
#           for letter in letters:
#               if letter not in group:
#                   break
#           else:
#               properties.append(name)
#               break
#   return properties
#
# catNames = getAlignProperties(alignment, AA_CATEGORIES)
# for i, category in enumerate(catNames):
#     print(i, category)