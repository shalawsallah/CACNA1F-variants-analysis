#!/usr/bin/env python3
# 10/05/17

import random

'''to calculate the randomness of the MODELLED mutants (through calculating change in residue-volume) in each set of variants placed in 4 equal bins according to the number of volume changes generated from differences in the 20 residues, ie, at 3 breaks: -42, 0, 42 '''

aa_volumes = {'G':66, 'A':92, 'S':99, 'C':106, 'T':122, 'D':125, 'P':129, 'N':135, 'V':142, 'E':155, 'Q':161, 'H':167, 'L':168, 'I':169, 'M':171, 'K':171, 'F':203, 'Y':203, 'R':225, 'W':240}
bin4_aa = ['I', 'F', 'Y', 'R', 'R', 'R', 'R', 'W'] # the mutated aa resulting in big loss in volume as a result of mutation
bin3_aa = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'N', 'N', 'S', 'I', 'E', 'P'] # the mutated aa resulting in small loss in volume as a result of mutation
bin2_aa = ['S', 'E', 'R', 'R', 'R', 'D', 'G', 'A', 'A', 'A', 'A'] # the mutated aa resulting in small gain in volume as a result of mutation
bin1_aa = ['L', 'L', 'L', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'S', 'P', 'D', 'C', 'C'] # the mutated aa resulting in big gain in volume as a result of mutation
all_aa = list(aa_volumes.keys())
# print(all_aa)
mutation_count = 54 # hgmd: highly conserved=25, moderately=26. Exac: highly conserved=33, moderately=53 (excluding 1 that has 0change in size).
observed_mutations_bin1 = 19 # bin1 size difference between mutated and mutant in mutations is < -42 Å3
observed_mutations_bin2 = 13 # bin 2 size difference between mutated and mutant in mutations is > -42 but < 0 Å3
observed_mutations_bin3 = 14 # bin 3 size difference between mutated and mutant in mutations is > 0 but < 42 Å3
observed_mutations_bin4 = 8 # bin 4 size difference between mutated and mutant in mutations is > 42 Å3

# bin1:
observed_mutations_in_bin1_fraction = float(observed_mutations_bin1) / float(mutation_count)
observed_mutations_outside_bin1_fraction = float(observed_mutations_bin2 + observed_mutations_bin3 + observed_mutations_bin4) / float(mutation_count)
print("Observed bin1 mutations fraction: ", observed_mutations_in_bin1_fraction)
print("Observed other bins mutations fraction: ", observed_mutations_outside_bin1_fraction)
# counter_higher: the number of times that the generated fraction, of "mutations with a large difference in volume", is greater than observed_fraction
# counter_lower: the number of times that the generated fraction of "mutations with a large difference in volume", is less than or equal to observed_fraction
counter_generated_mutations_outside_bin1_higher = 0
counter_generated_mutations_outside_bin1_lower = 0
counter_generated_mutations_in_bin1_higher = 0
counter_generated_mutations_in_bin1_lower = 0
n_iterations = 10000
for j in range(n_iterations):
    randomly_generated_aa = []# Repetitions allowed
    for m in range(len(all_aa)):
        randomly_generated_aa.append(all_aa[random.randint(0, len(all_aa) - 1)]) # random.randint between 1 and 20 (length of all_aa). Then, all_aa[]: random integers that correspond to elements (assigned randomly) that are within all_aa[]
        # print('randomly generated aa: ', randomly_generated_aa)
    in_bin1_counter = 0 # number of randomly generated mutations in bin1
    outside_bin1_counter = 0 # number of randomly generated mutations outside bin1
    cutoff_mutated_mutant_size_difference = -42 # separating the bins

    for amino in range(len(bin1_aa)):
        for each_aa in range(len(randomly_generated_aa)):
            volume_difference = aa_volumes[bin1_aa[amino]] - aa_volumes[randomly_generated_aa[each_aa]] 
            # print(volume_difference)
            if volume_difference < cutoff_mutated_mutant_size_difference:
                in_bin1_counter += 1
            elif volume_difference >= cutoff_mutated_mutant_size_difference:
                outside_bin1_counter += 1
            else:
                print('problem')

    generated_mutations_outside_bin1_fraction = float(outside_bin1_counter) / float(outside_bin1_counter + in_bin1_counter)
    generated_mutations_in_bin1_fraction = float(in_bin1_counter) / float(outside_bin1_counter + in_bin1_counter)

    if generated_mutations_outside_bin1_fraction > observed_mutations_outside_bin1_fraction:
        counter_generated_mutations_outside_bin1_higher += 1
    elif generated_mutations_outside_bin1_fraction <= observed_mutations_outside_bin1_fraction:
        counter_generated_mutations_outside_bin1_lower += 1
    else:
       print("problem1")
    if generated_mutations_in_bin1_fraction > observed_mutations_in_bin1_fraction:
        counter_generated_mutations_in_bin1_higher += 1
    elif generated_mutations_in_bin1_fraction <= observed_mutations_in_bin1_fraction:
        counter_generated_mutations_in_bin1_lower += 1
    else:
       print("problem2")
print('number of randomly generated aa: ', len(randomly_generated_aa))
print('both counters: ', in_bin1_counter + outside_bin1_counter)

# print ("counter_generated_mutations_in_bin1_lower:") # this counter = counter_generated_mutations_outside_bin1_higher
# print (counter_generated_mutations_in_bin1_lower)
# print ("counter_generated_mutations_in_bin1_higher:")# this counter = counter_generated_mutations_outside_bin1_lower
# print (counter_generated_mutations_in_bin1_higher)
# print ("counter_generated_mutations_outside_bin1_lower:")
# print (counter_generated_mutations_outside_bin1_lower)
# print ("counter_generated_mutations_outside_bin1_higher:")
# print (counter_generated_mutations_outside_bin1_higher)

# In regard to possible mutation outcomes:
expected_mutations_in_bin1 = observed_mutations_in_bin1_fraction * (observed_mutations_bin2 + observed_mutations_bin3 + observed_mutations_bin4) / mutation_count
print('Expected mutations in bin1: ', expected_mutations_in_bin1)

expected_mutations_outside_bin1 = observed_mutations_outside_bin1_fraction * (observed_mutations_bin1) / mutation_count
print('Expected mutations outside bin1: ', expected_mutations_outside_bin1)

p_mutations_outside_bin1 = float(counter_generated_mutations_outside_bin1_higher) / float(counter_generated_mutations_outside_bin1_lower + counter_generated_mutations_outside_bin1_higher)
p_mutations_in_bin1 = float(counter_generated_mutations_in_bin1_higher) / float(counter_generated_mutations_in_bin1_lower + counter_generated_mutations_in_bin1_higher) #opposite of above

print("outside bin1 mutations p-value: ", p_mutations_outside_bin1)
print("in bin1 mutations p-value:" , p_mutations_in_bin1)
print('#######################################')

# bin2:
observed_mutations_in_bin2_fraction = float(observed_mutations_bin2) / float(mutation_count)
observed_mutations_outside_bin2_fraction = float(observed_mutations_bin1 + observed_mutations_bin3 + observed_mutations_bin4) / float(mutation_count)

# counter_higher: the number of times that the generated fraction, of "mutations with a large difference in volume", is greater than observed_fraction
# counter_lower: the number of times that the generated fraction of "mutations with a large difference in volume", is less than or equal to observed_fraction
counter_generated_mutations_outside_bin2_higher = 0
counter_generated_mutations_outside_bin2_lower = 0

counter_generated_mutations_in_bin2_higher = 0
counter_generated_mutations_in_bin2_lower = 0

print("Observed bin2 mutations fraction:")
print(observed_mutations_in_bin2_fraction)

print("Observed other bins mutations fraction:")
print(observed_mutations_outside_bin2_fraction)

for j in range(n_iterations):

    randomly_generated_aa = []
    # Repetitions allowed
    for m in range(len(all_aa)):
        randomly_generated_aa.append(all_aa[random.randint(0, len(all_aa) - 1)]) # random.randint between 1 and 20 (length of all_aa). Then, all_aa[]: random integers that correspond to elements (assigned randomly) that are within all_aa[]
        # print(randomly_generated_aa)
    # print(randomly_generated_aa[0])
    # print(aa_volumes[randomly_generated_aa[0]])

    # in_bin_counter : number of randomly generated mutations in bin
    # outside_bin_counter: number of randomly generated mutations outside bin
    in_bin2_counter = 0
    outside_bin2_counter = 0
    cutoff_mutated_mutant_size_difference = -42 # separating the bins

    for amino in range(len(bin2_aa)):
        for each_aa in range(len(randomly_generated_aa)):
            volume_difference = aa_volumes[bin2_aa[amino]] - aa_volumes[randomly_generated_aa[each_aa]] # aa_volumes[] calls the values of each amino acid in bin2 and those generated before subtraction
            # print(volume_difference)
            if volume_difference < 0 and volume_difference >= cutoff_mutated_mutant_size_difference:
                in_bin2_counter += 1
            elif volume_difference < cutoff_mutated_mutant_size_difference or volume_difference >= 0:
                outside_bin2_counter += 1
            else:
                print('problem')

    generated_mutations_outside_bin2_fraction = float(outside_bin2_counter) / float(outside_bin2_counter + in_bin2_counter)
    generated_mutations_in_bin2_fraction = float(in_bin2_counter) / float(outside_bin2_counter + in_bin2_counter)

    # print('generated_mutations_outside_bin4_fraction: ')
    # print (generated_mutations_outside_bin4_fraction)
    # print('generated_mutations_in_bin4_fraction: ')
    # print (generated_mutations_in_bin4_fraction)

    if generated_mutations_outside_bin2_fraction > observed_mutations_outside_bin2_fraction:
        counter_generated_mutations_outside_bin2_higher += 1
    elif generated_mutations_outside_bin2_fraction <= observed_mutations_outside_bin2_fraction:
        counter_generated_mutations_outside_bin2_lower += 1
    else:
       print("problem1")

    if generated_mutations_in_bin2_fraction > observed_mutations_in_bin2_fraction:
        counter_generated_mutations_in_bin2_higher += 1
    elif generated_mutations_in_bin2_fraction <= observed_mutations_in_bin2_fraction:
        counter_generated_mutations_in_bin2_lower += 1
    else:
       print("problem2")

print('number of randomly generated aa: ')
print(len(randomly_generated_aa))
print('both counters: ')
print (in_bin2_counter + outside_bin2_counter)

# print ("counter_generated_mutations_in_bin2_lower:") # this counter = counter_generated_mutations_outside_bin2_higher
# print (counter_generated_mutations_in_bin2_lower)
# print ("counter_generated_mutations_in_bin2_higher:")# this counter = counter_generated_mutations_outside_bin2_lower
# print (counter_generated_mutations_in_bin2_higher)
#
# print ("counter_generated_mutations_outside_bin2_lower:")
# print (counter_generated_mutations_outside_bin2_lower)
# print ("counter_generated_mutations_outside_bin2_higher:")
# print (counter_generated_mutations_outside_bin2_higher)

# In regard to possible mutation outcomes:
expected_mutations_in_bin2 = observed_mutations_in_bin2_fraction * (observed_mutations_bin1 + observed_mutations_bin3 + observed_mutations_bin4) / mutation_count
print('Expected mutations in bin2: ')
print(expected_mutations_in_bin2)

expected_mutations_outside_bin2 = observed_mutations_outside_bin2_fraction * (observed_mutations_bin2) / mutation_count
print('Expected mutations outside bin2: ')
print(expected_mutations_outside_bin2)

p_mutations_outside_bin2 = float(counter_generated_mutations_outside_bin2_higher) / float(counter_generated_mutations_outside_bin2_lower + counter_generated_mutations_outside_bin2_higher)
p_mutations_in_bin2 = float(counter_generated_mutations_in_bin2_higher) / float(counter_generated_mutations_in_bin2_lower + counter_generated_mutations_in_bin2_higher)

print("outside bin2 mutations p-value:")
print(p_mutations_outside_bin2)
print("in bin2 mutations p-value:")
print(p_mutations_in_bin2)
print('#######################################')

# bin3:
observed_mutations_in_bin3_fraction = float(observed_mutations_bin3) / float(mutation_count)
observed_mutations_outside_bin3_fraction = float(observed_mutations_bin1 + observed_mutations_bin2 + observed_mutations_bin4) / float(mutation_count)

# counter_higher: the number of times that the generated fraction, of "mutations with a large difference in volume", is greater than observed_fraction
# counter_lower: the number of times that the generated fraction of "mutations with a large difference in volume", is less than or equal to observed_fraction
counter_generated_mutations_outside_bin3_higher = 0
counter_generated_mutations_outside_bin3_lower = 0

counter_generated_mutations_in_bin3_higher = 0
counter_generated_mutations_in_bin3_lower = 0

print("Observed bin3 mutations fraction:")
print(observed_mutations_in_bin3_fraction)

print("Observed other bins mutations fraction:")
print(observed_mutations_outside_bin3_fraction)

for j in range(n_iterations):

    randomly_generated_aa = []
    # Repetitions allowed
    for m in range(len(all_aa)):
        randomly_generated_aa.append(all_aa[random.randint(0, len(all_aa) - 1)]) # random.randint between 1 and 20 (length of all_aa). Then, all_aa[]: random integers that correspond to elements (assigned randomly) that are within all_aa[]
        # print(randomly_generated_aa)
    # print(randomly_generated_aa[0])
    # print(aa_volumes[randomly_generated_aa[0]])

    # in_bin_counter : number of randomly generated mutations in bin
    # outside_bin_counter: number of randomly generated mutations outside bin
    in_bin3_counter = 0
    outside_bin3_counter = 0
    cutoff_mutated_mutant_size_difference = 42 # separating the bins

    for amino in range(len(bin3_aa)):
        for each_aa in range(len(randomly_generated_aa)):
            volume_difference = aa_volumes[bin3_aa[amino]] - aa_volumes[randomly_generated_aa[each_aa]] # aa_volumes[] calls the values of each amino acid in bin2 and those generated before subtraction
            # print(volume_difference)
            if volume_difference > 0 and volume_difference <= cutoff_mutated_mutant_size_difference:
                in_bin3_counter += 1
            elif volume_difference > cutoff_mutated_mutant_size_difference or volume_difference <= 0:
                outside_bin3_counter += 1
            else:
                print('problem')

    generated_mutations_outside_bin3_fraction = float(outside_bin3_counter) / float(outside_bin3_counter + in_bin3_counter)
    generated_mutations_in_bin3_fraction = float(in_bin3_counter) / float(outside_bin3_counter + in_bin3_counter)

    # print('generated_mutations_outside_bin4_fraction: ')
    # print (generated_mutations_outside_bin4_fraction)
    # print('generated_mutations_in_bin4_fraction: ')
    # print (generated_mutations_in_bin4_fraction)

    if generated_mutations_outside_bin3_fraction > observed_mutations_outside_bin3_fraction:
        counter_generated_mutations_outside_bin3_higher += 1
    elif generated_mutations_outside_bin3_fraction <= observed_mutations_outside_bin3_fraction:
        counter_generated_mutations_outside_bin3_lower += 1
    else:
       print("problem1")

    if generated_mutations_in_bin3_fraction > observed_mutations_in_bin3_fraction:
        counter_generated_mutations_in_bin3_higher += 1
    elif generated_mutations_in_bin3_fraction <= observed_mutations_in_bin3_fraction:
        counter_generated_mutations_in_bin3_lower += 1
    else:
       print("problem2")

print('number of randomly generated aa: ')
print(len(randomly_generated_aa))

print('both counters: ')
print (in_bin3_counter + outside_bin3_counter)

# print ("counter_generated_mutations_in_bin3_lower:") # this counter = counter_generated_mutations_outside_bin3_higher
# print (counter_generated_mutations_in_bin3_lower)
# print ("counter_generated_mutations_in_bin3_higher:")# this counter = counter_generated_mutations_outside_bin3_lower
# print (counter_generated_mutations_in_bin3_higher)
#
# print ("counter_generated_mutations_outside_bin3_lower:")
# print (counter_generated_mutations_outside_bin3_lower)
# print ("counter_generated_mutations_outside_bin3_higher:")
# print (counter_generated_mutations_outside_bin3_higher)

# In regard to possible mutation outcomes:
expected_mutations_in_bin3 = observed_mutations_in_bin3_fraction * (observed_mutations_bin1 + observed_mutations_bin2 + observed_mutations_bin4) / mutation_count
print('Expected mutations in bin3: ')
print(expected_mutations_in_bin3)

expected_mutations_outside_bin3 = observed_mutations_outside_bin3_fraction * (observed_mutations_bin3) / mutation_count
print('Expected mutations outside bin3: ')
print(expected_mutations_outside_bin3)

p_mutations_outside_bin3 = float(counter_generated_mutations_outside_bin3_higher) / float(counter_generated_mutations_outside_bin3_lower + counter_generated_mutations_outside_bin3_higher)
p_mutations_in_bin3 = float(counter_generated_mutations_in_bin3_higher) / float(counter_generated_mutations_in_bin3_lower + counter_generated_mutations_in_bin3_higher)

print("outside bin3 mutations p-value:")
print(p_mutations_outside_bin3)

print("in bin3 mutations p-value:")
print(p_mutations_in_bin3)

print('#######################################')

# bin4:
observed_mutations_in_bin4_fraction = float(observed_mutations_bin4) / float(mutation_count)
observed_mutations_outside_bin4_fraction = float(observed_mutations_bin1 + observed_mutations_bin2 + observed_mutations_bin3) / float(mutation_count)

# counter_higher: the number of times that the generated fraction, of "mutations with a large difference in volume", is greater than observed_fraction
# counter_lower: the number of times that the generated fraction of "mutations with a large difference in volume", is less than or equal to observed_fraction
counter_generated_mutations_outside_bin4_higher = 0
counter_generated_mutations_outside_bin4_lower = 0

counter_generated_mutations_in_bin4_higher = 0
counter_generated_mutations_in_bin4_lower = 0

print("Observed bin4 mutations fraction:")
print(observed_mutations_in_bin4_fraction)

print("Observed other bins mutations fraction:")
print(observed_mutations_outside_bin4_fraction)

for j in range(n_iterations):

    randomly_generated_aa = []
    # Repetitions allowed
    for m in range(len(all_aa)):
        randomly_generated_aa.append(all_aa[random.randint(0, len(all_aa) - 1)]) # random.randint between 1 and 20 (length of all_aa). Then, all_aa[]: random integers that correspond to elements (assigned randomly) that are within all_aa[]
        # print(randomly_generated_aa)
    # print(randomly_generated_aa[0])
    # print(aa_volumes[randomly_generated_aa[0]])

    # in_bin_counter : number of randomly generated mutations in bin
    # outside_bin_counter: number of randomly generated mutations outside bin
    in_bin4_counter = 0
    outside_bin4_counter = 0
    cutoff_mutated_mutant_size_difference = 42 # separating the bins

    for amino in range(len(bin4_aa)):
        for each_aa in range(len(randomly_generated_aa)):
            volume_difference = aa_volumes[bin4_aa[amino]] - aa_volumes[randomly_generated_aa[each_aa]] # aa_volumes[] calls the values of each amino acid in bin2 and those generated before subtraction
            # print(volume_difference)
            if volume_difference > cutoff_mutated_mutant_size_difference:
                in_bin4_counter += 1
            elif volume_difference <= cutoff_mutated_mutant_size_difference:
                outside_bin4_counter += 1
            else:
                print('problem')

    generated_mutations_outside_bin4_fraction = float(outside_bin4_counter) / float(outside_bin4_counter + in_bin4_counter)
    generated_mutations_in_bin4_fraction = float(in_bin4_counter) / float(outside_bin4_counter + in_bin4_counter)

    # print('generated_mutations_outside_bin4_fraction: ')
    # print (generated_mutations_outside_bin4_fraction)
    # print('generated_mutations_in_bin4_fraction: ')
    # print (generated_mutations_in_bin4_fraction)

    if generated_mutations_outside_bin4_fraction > observed_mutations_outside_bin4_fraction:
        counter_generated_mutations_outside_bin4_higher += 1
    elif generated_mutations_outside_bin4_fraction <= observed_mutations_outside_bin4_fraction:
        counter_generated_mutations_outside_bin4_lower += 1
    else:
       print("problem1")

    if generated_mutations_in_bin4_fraction > observed_mutations_in_bin4_fraction:
        counter_generated_mutations_in_bin4_higher += 1
    elif generated_mutations_in_bin4_fraction <= observed_mutations_in_bin4_fraction:
        counter_generated_mutations_in_bin4_lower += 1
    else:
       print("problem2")

print('number of randomly generated aa: ')
print(len(randomly_generated_aa))

print('both counters: ')
print (in_bin4_counter + outside_bin4_counter)

# print ("counter_generated_mutations_in_bin4_lower:") # this counter = counter_generated_mutations_outside_bin4_higher
# print (counter_generated_mutations_in_bin4_lower)
# print ("counter_generated_mutations_in_bin4_higher:")# this counter = counter_generated_mutations_outside_bin4_lower
# print (counter_generated_mutations_in_bin4_higher)
#
# print ("counter_generated_mutations_outside_bin4_lower:")
# print (counter_generated_mutations_outside_bin4_lower)
# print ("counter_generated_mutations_outside_bin4_higher:")
# print (counter_generated_mutations_outside_bin4_higher)

# In regard to possible mutation outcomes:
expected_mutations_in_bin4 = observed_mutations_in_bin4_fraction * (observed_mutations_bin1 + observed_mutations_bin2 + observed_mutations_bin3) / mutation_count
print('Expected mutations in bin4: ')
print(expected_mutations_in_bin4)

expected_mutations_outside_bin4 = observed_mutations_outside_bin4_fraction * (observed_mutations_bin4) / mutation_count
print('Expected mutations outside bin4: ')
print(expected_mutations_outside_bin4)

p_mutations_outside_bin4 = float(counter_generated_mutations_outside_bin4_higher) / float(counter_generated_mutations_outside_bin4_lower + counter_generated_mutations_outside_bin4_higher)
p_mutations_in_bin4 = float(counter_generated_mutations_in_bin4_higher) / float(counter_generated_mutations_in_bin4_lower + counter_generated_mutations_in_bin4_higher)

print("outside bin4 mutations p-value:")
print(p_mutations_outside_bin4)

print("in bin4 mutations p-value:")
print(p_mutations_in_bin4)