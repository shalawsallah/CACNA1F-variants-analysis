#! /usr/bin/env python3
import pandas as pd
import os
import re
'''To batch-query the database using Genomic coordinates obtained from VariantValidator.
Before running the script SORT the spreadsheet (A to Z for forward strands & reverse for genes on minus strands) based on coordinates, variants, and nucleotide change to make the alignment with the corresponding dbNSFP predictions later easier. this alignment is needed to carry over the class pathogenicicty of the corresponding variant from the analysis sheet, & because these  predictions automatically exclude duplications'''

gene = 'cryba1', 'cryba2', 'cryba4', 'crybb1', 'crybb2', 'crybb3', 'crygc', 'crygd', 'crygs'
chromosome = '17', '2', '22', '22', '22', '22', '2', '2', '3'
# strand = '1' # from Ensembl:1 if forward or -1 if reverse
#######################################################
f_vars = pd.read_excel('/Users/mdefsss2/other_genes/crystallins/variant_validator_output.xlsx')
output_file = open('/Users/mdefsss2/dbNSFP4/crystallins.in', "w+")
file_contains_gene_NM = pd.read_excel('/Users/mdefsss2/other_genes/genes_list.xlsx')
coordinates = f_vars['GRCh38_POS']
ref_allele = f_vars['GRCh38_REF']
alt_allele = f_vars['GRCh38_ALT']
variants = f_vars['Input']
genes = f_vars['Gene_Symbol']
for g in range(len(gene)):
    refSeq = ''
    for i in range(len(file_contains_gene_NM)):
        if file_contains_gene_NM['gene'][i] == gene[g]:
            refSeq += file_contains_gene_NM['transcript'][i]

    '''for dbNSFP input'''
    for i in range(len(variants)):
        if genes[i] == gene[g].upper():
            output_file.write(str(chromosome[g])+' '+str(coordinates[i])+' '+str(ref_allele[i])+' '+str(alt_allele[i])+'\n')
    # print(str(chromosome)+' '+str(coordinates[i])+' '+str(ref_allele[i])+' '+str(alt_allele[i]))
output_file.close()
print('results written to file \n variants submitted to dbNSFP..')
exit(True)
# input ='java search_dbNSFP40a -i '+gene+'.in -o '+gene+'.out -w 2-6,9,12-13,15,37,46,49,53,58,61,64,67,69,72,76,79,81,86,88,90,93,104-105,107,120,137,139,145,151,226,229,365-367,373,394-395,412,415,417,431-433,438,441,450-452,455-456,459,464'
# input ='java search_dbNSFP40a -i '+gene[g]+'.in -o '+gene[g]+'.out'
input ='java search_dbNSFP40a -i crystallins.in -o crystallins.out'
print(input)
os.chdir("/Users/mdefsss2/dbNSFP4/")
os.system(input)