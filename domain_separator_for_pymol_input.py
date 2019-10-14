# edited on 07/05/17
# pdb file used excludes the CTD as well. it only includes domains 1-4

# created 25/04/17
# producing the format needed by Pymol of the protein residues to be separated into spatial domains (top & bottom):
# Working on the following script (which is angleCalculator_lowerUpperDomains.py):

# ! /usr/bin/env python3
# 12/04/17
# Reading xyz coordinates from a pdb file to measure angles between Calpha atoms,
# in order to compute spatial distribution of mutations in cacna1f:

import math
import sys
from Bio.PDB import *
import numpy as np
import cmd

parser = PDBParser()
# Parse the structure into a PDB.Structure object
pdb_code = open('/Users/mdefsss2/rpgr/4jhn_reduced.pdb')  # file exludes CTD and N_terminal
pdb_path = "/Users/mdefsss2/rpgr/4jhn_reduced.pdb"
struct = parser.get_structure(pdb_code, pdb_path)
pdb = open('/Users/mdefsss2/rpgr/4jhn_reduced.pdb')
residue_coordinates = []
res_one = []
x_coord = []
y_coord = []
z_coord = []
counter = 0

for line in pdb:
	##    if line.startswith('ATOM'):
	if line[:6] == 'ATOM  ' and line[12:16] == ' CA ':
		##        print(line)
		atom_number = int(line[6:11].strip())
		atom_name = line[12:16].strip()
		residue_name = line[17:20].strip()
		residue_number = int(line[22:26].strip())
		# print (residue_number)
		x_coord.append(float(line[30:38].strip()))
		y_coord.append(float(line[38:46].strip()))
		z_coord.append(float(line[46:54].strip()))
		counter += 1
		residue_coordinates.append([atom_name, residue_name, residue_number, x_coord, y_coord, z_coord])
		# print(len(residue_coordinates))

# may it be necessary to specify a model, and a chain? in which case the following maybe needed: model = structure[0]

# Construct list of vectors
list_of_vectors = []
# list_of_residues = []
for model in struct:
	# print(model)
	for chain in model:
		# print(type(chain))
		for residue in chain:
			# print(residue.get_id()[1])
			# print(type(residue))
			for atom in residue:
				# print(atom)
				id = atom.id
				# print(id)

				# first way, including a counter which can be included along the list_of_vectors
				#  this can be used to refer to the number the 2 constant atoms used to calculate angles
				#                 if id == 'CA':
				#                     counter += 1
				#                     vector = atom.get_vector()
				#                     list_of_vectors.append([counter, vector])
				#                     # print(type(list_of_vectors))
				#                     ## enumerate() could be used?
				#                     # x = list(enumerate(list_of_vectors))
				#                     # print(x)
				#                     # print(id)
				# and to iterate would be:
				# for i in range(len(list_of_vectors)):
				#     print(list_of_vectors[i][0], ": ", i)

				# the other way, excluding a counter:
				if id == 'CA':
					vector = atom.get_vector()
					# print(vector)
					list_of_vectors.append(
						[vector, residue.get_id()[1]])  # residue.get_id() contains more than one element, hence [1]
					#                     list_of_residues.append(residue.get_id()[1]) # another way of adding residue.id to a list
					# print(list_of_vectors)
					# x = list(enumerate(list_of_vectors))

# calculating centre-mass of the protein to use in angle calculation
mean_x_coord = sum(x_coord) / len(x_coord)
# print(mean_x_coord)
mean_y_coord = sum(y_coord) / len(y_coord)
# print(mean_y_coord)
mean_z_coord = sum(z_coord) / len(z_coord)
# print(mean_z_coord)
centre = [mean_x_coord, mean_y_coord, mean_z_coord]
# print(centre)
# print(type(centre))


pi_by_two = np.pi / 2.0  # 90 degrees from radian
# Use the vector representation of the atomic coordinates, and the calc angle function from the Vector module:
vector1 = list_of_vectors[304][0]  # residue number 383 assigned based on visualisation. it is situated in the middle, below the pore
# However, since the file starts from residue number 79 (domain 1), 383 does not correspond to the residue of interest (it would be 383 minus 78, and accounting for [0] element)
print(vector1)
vector2 = Vector(centre)  # the centre-mass is assigned as the centre of the protein
print(vector2)

# Iterate through the list to calculate angles between vectors 1, 2 and every other ca atom as vector3:
counter_bottom_half = 0
counter_top_half = 0

string_of_bottom_residues = ''
first_bottom_residue = 0
last_bottom_residue = 0
bottom_initial = True

string_of_top_residues = ''
first_top_residue = 0
last_top_residue = 0
top_initial = True

# the first residue in the list is the same:
first_top_residue = list_of_vectors[1]
last_top_residue = list_of_vectors[1]

for res in range(len(list_of_vectors)):
	# print(res, ": ", list_of_vectors[res])
	vector3 = list_of_vectors[res][0]
	angle = calc_angle(vector1, vector2, vector3)
	# print(angle)
	# print(vector3)

	if angle <= pi_by_two:
		counter_bottom_half += 1
		# string_of_bottom_residues += str(list_of_vectors[res][1])
		#list_of_bottom_residues.append(list_of_vectors[res][1])

        # Bottom part:
		if bottom_initial:  # First iteration, for the first residue in the list
			first_bottom_residue = list_of_vectors[res][1]
			last_bottom_residue = list_of_vectors[res][1]
			bottom_initial = False
		else:
			if list_of_vectors[res][1] == last_bottom_residue + 1:
				last_bottom_residue = list_of_vectors[res][1]
			else:
				if first_bottom_residue != last_bottom_residue:
					string_of_bottom_residues = string_of_bottom_residues + str(first_bottom_residue) + "-" + str(
						last_bottom_residue) + '+' # adding a consecutive numbers (in a range) to the string of residues
				else:
					string_of_bottom_residues = string_of_bottom_residues + str(first_bottom_residue) + '+'
				first_bottom_residue = list_of_vectors[res][1]
				last_bottom_residue = list_of_vectors[res][1]
	else:
		counter_top_half += 1
		# string_of_top_residues += '+' + str(list_of_vectors[res][1])
        # Top part:
		if top_initial:  # First iteration, for the first residue in the list
			first_top_residue = list_of_vectors[res][1]
			last_top_residue = list_of_vectors[res][1]
			top_initial = False
		else:
			if list_of_vectors[res][1] == last_top_residue + 1:
				last_top_residue = list_of_vectors[res][1]
			else:
				if first_top_residue != last_top_residue:
					string_of_top_residues = string_of_top_residues + str(first_top_residue) + "-" + str(
						last_top_residue) + '+' # adding a consecutive numbers (in a range) to the string of residues
				else:
					string_of_top_residues = string_of_top_residues + str(first_top_residue) + '+'
				first_top_residue = list_of_vectors[res][1]
				last_top_residue = list_of_vectors[res][1]
# print(list_of_bottom_residues)

# for last set of consecutive numbers of on domain:
if first_bottom_residue != last_bottom_residue:
	string_of_bottom_residues = string_of_bottom_residues + str(first_bottom_residue) + "-" + str(last_bottom_residue)
else:
	string_of_bottom_residues = string_of_bottom_residues + str(first_bottom_residue)
print('bottom residues: ', string_of_bottom_residues)

# for last set of consecutive numbers of the other domain:
if first_top_residue != last_top_residue:
	string_of_top_residues = string_of_top_residues + str(first_top_residue) + "-" + str(last_top_residue)
else:
	string_of_top_residues = string_of_top_residues + str(first_top_residue)
print('top residues: ', string_of_top_residues)

# first_top_residue = string_of_top_residues[0]
# last_top_residue = string_of_top_residues[0]
#
# for i in range(len(string_of_top_residues)):
#
# 	if i == 0:  # First number
# 		first_top_residue = string_of_top_residues[i]
# 		last_top_residue = string_of_top_residues[i]
# 	else:
# 		if string_of_top_residues[i] == str(last_top_residue + 1):
# 			last_top_residue = string_of_top_residues[i]
# 		else:
# 			if first_top_residue != last_top_residue:
# 				string_pf_top_residues = string_of_top_residues + str(first_top_residue) + "-" + str(last_top_residue) + '+'
# 			else:
# 				string_of_top_residues = string_of_top_residues + str(first_top_residue) + '+'
# 				first_top_residue = string_of_top_residues[i]
# 				last_top_residue = string_of_top_residues[i]
#
# if first_top_residue != last_top_residue:
# 	string_of_top_residues = string_of_top_residues + str(first_top_residue) + "-" + str(last_top_residue)
# else:
# 	string_of_top_residues = str(string_of_top_residues) + str(first_top_residue)
# print(string_of_top_residues)


print('bottom half Calphas: ' + str(counter_bottom_half))
print('top half Calphas: ' + str(counter_top_half))
# print('bottom residues: ')
# print(string_of_bottom_residues[1:])
# print('top residues: ')
# print(string_of_top_residues[1:])
