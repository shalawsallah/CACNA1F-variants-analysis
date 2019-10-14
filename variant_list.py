#!/usr/bin/python

'''28/10/18

PYMOL SCRIPT 2.

inputing the variant numbers to output in the format required in PyMol command line.
beforehand, residue numbers are extracted using 'extract_resi_numbers.py' & copied into 'pymol_test_residues.txt' for each gene separately
'''

f1=open('/Users/mdefsss2/PycharmProjects/spatial_distribution/pymol_visualisation/pathogenic_CACNA1F_resNumbers.txt')
array1 = ''
filecontents1 = f1.read().split()
for i in filecontents1:
            x= (i+ '+')
            array1 = array1 + x

print(array1)
f1.close()

print('################################')

f2=open('/Users/mdefsss2/PycharmProjects/spatial_distribution/pymol_visualisation/non_pathogenic_CACNA1F_resNumbers.txt')
array2 = ''
filecontents2 = f2.read().split()
for i in filecontents2:
            x= (i+ '+')
            array2 = array2 + x

print(array2)
f2.close()