#!/usr/bin/env python3
# 07/05/17

import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
'''To calculate side chain differences between mutated and mutants into a matrix of amino acid & a plot a histogram of its distribution'''
# aa = [Gly,Ala,Ser,Cys, Thr, Asp, Pro, Asn, Val, Glu, Gln, His, Leu, Ile, Lys, Met, Tyr, Phe, Arg, Trp]
aa_1 = [66, 92, 99, 106, 122, 125, 129, 135, 142, 155, 161, 167, 168, 169, 171, 171, 203, 203, 225, 240]
aa_2 = [66, 92, 99, 106, 122, 125, 129, 135, 142, 155, 161, 167, 168, 169, 171, 171, 203, 203, 225, 240]
string = ''
vol_difference_20x20 = []
for i in range (len(aa_1)):
    for j in range (len(aa_2)):
        '''in matrix format'''
#         string += str(aa_1[i] - aa_2[j]) + ' '
#     string = string[:-1] + '\n'
# print (string)

        vol_difference_20x20.append( aa_1[i] - aa_2[j] )
print(vol_difference_20x20) # a list of the differences

'''the following places the numbers in the above list into four bins of about 100 residues, including those with zero difference'''
in_bin1_counter = 0
in_bin2_counter = 0
in_bin3_counter = 0
in_bin4_counter = 0
zero_counter = 0
bin1_differences = []
for val in vol_difference_20x20:
    if val < -42: # when having the ranges (bins in the case) in order, there is no need for a range, ie. val < -42 and val > 0
        in_bin1_counter += 1
        bin1_differences.append(val)
    elif val < 0:
        in_bin2_counter += 1
    elif val == 0:
        zero_counter += 1
    elif val <= 42:
        in_bin3_counter += 1
    elif val >42:
        in_bin4_counter += 1
print('in_bin1_counter ' + str(in_bin1_counter))
print('in_bin2_counter ' + str(in_bin2_counter))
print('in_bin3_counter ' + str(in_bin3_counter))
print('in_bin4_counter ' + str(in_bin4_counter))
print('zero_counter ' + str(zero_counter))

## Calculating the fraction of variants in each bin, based on difference in residue volume caused by a mutation
in_bin1_res_vol_difference_fraction = in_bin1_counter / len(vol_difference_20x20)
print('number of bin1 variants ' + str(in_bin1_res_vol_difference_fraction))
in_bin2_res_vol_difference_fraction = in_bin2_counter / len(vol_difference_20x20)
print('number of bin2 variants ' + str(in_bin2_res_vol_difference_fraction))
in_bin3_res_vol_difference_fraction = in_bin3_counter / len(vol_difference_20x20)
print('number of bin3 variants ' + str(in_bin3_res_vol_difference_fraction))
in_bin4_res_vol_difference_fraction = in_bin4_counter / len(vol_difference_20x20)
print('number of bin4 variants ' + str(in_bin4_res_vol_difference_fraction))

# generating a histogram for the distribution of the resulting difference in the volumes of residues, due to mutation:
bins = 17 # 0r 15
plt.hist(vol_difference_20x20, bins, rwidth=0.98)  # plot

plt.xlabel('20*20 side chain differences')
plt.ylabel('frequency')
plt.title('side chain-volume change distribution')
# plt.legend()
plt.show()

# alpha = 0.05 #( this is 95% of confidence)
'''Test statistical significance'''
bin1_pvalue = stats.mannwhitneyu( vol_difference_20x20, bin1_differences, alternative='two-sided' )
print(bin1_pvalue)
# result = "\n==> Not Statistically Significant (p-value: %f)" % pvalue
#if pvalue < alpha:
    # result = "\n==> Statistically Significant (p-value: %f)" % pvalue
 #   print ('significant')
# print(result) # the second element of 'result' is the p value. ie. result = stats.mannwhitneyu( aa_volumes, bin1_aa, alternative='two-sided' )[1]

sns.distplot(vol_difference_20x20)
plt.show()