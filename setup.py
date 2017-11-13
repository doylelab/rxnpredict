from scipy.stats.stats import pearsonr
from subprocess import call
import numpy as np
import math
import sys 
import glob 
import subprocess 
import os 
import shutil 
import collections
import time
import platform
import csv

# USER: change to path of Spartan executable file (use \\ and exclude .exe)
spartan_path = 'C:\\Program Files\\Wavefunction\\Spartan14v114'

localtime = time.asctime(time.localtime(time.time()))
print('*' * 100)
print('')
print('Time run:', localtime) 
print('')

header = [
	"C OPT B3LYP 6-31G* FINDBEST=MMFF FREQ NMR POLAR\n",
	"C CONVERGE SCFTOLERANCE=VERYHIGH GRADIENTTOLERANCE=0.000005\n",
	"DISTANCETOLERANCE=0.00002 BIGGRID\n"
]
footer = [
	"NOCONFORMER\n",
	"BEGINPROPIN\n",
	"C DOQSAR PRINTMO THERMO PRINTFREQ PRINTCHG\n", 
	"C QSAR NBO DIPOLE POP PROPPRINTLEV=2 PROPRERUN\n",
	" PRINTVIBCOORDS PROP:IGNORE_WARN\n",
	"ENDPROPIN\n",
	"BEGINPREFERENCES\n", 
	" MM:CONF_SELECTION_RULE=2\n",
	"ENDPREFERENCES\n"
]
label = "*"
atomic_weights = {
	'H': 1.0079, 
	'C': 12.0107, 
	'N': 14.0067, 
	'O': 15.9994, 
	'F': 18.9984,
	'P': 30.9737,
	'S': 32.065, 
	'Cl': 35.453, 
	'Br': 79.904,
	'I': 126.9045,
}


def split(text):
	"""Takes a list containing lines of text and splits each line to 
	make a list of lists.
	"""
	split_text = []
	for line in text: 
		split_text.append(line.split()) 
	return split_text

def in_between(text_list, start_keyword, end_keyword): 
	"""Takes a list of text and returns the list containing entries between
	start_keyword and end_keyword (inclusive).
	"""
	output = []
	parsing = False
	for line in text_list:
		if start_keyword in line: 
			parsing = True 
		if parsing:
			output.append(line) 
		if end_keyword in line: 
			parsing = False
	return output

def extract_table(text_list, shared, vname, key_col, value_col):
	"""Takes a list of text and subsets to the rows that contain the label.
	Then returns a dictionary with keys from key_col and values from value_col.
	Used for extracting multiple values from table in output file.

	Needs to take a list of shared keys and subset to those
	""" 
	keys = [] 
	values = []
	for line in text_list:
		if any(s in shared for s in line):
			keys.append(line[key_col-1])
			values.append(line[value_col-1])
	labeled_keys = [s + "_" + vname for s in keys]
	dictionary = dict(zip(labeled_keys, values))
	return dictionary

def extract_line(text_list, keyword, vname, col):
	"""Takes a list of text and for every line containing the keyword,
	creates a dictionary where vname is key and the word at col is the value.
	Used for extracting one value from output file.
	"""
	dictionary = {}
	for line in text_list:
		if keyword in line:
			dictionary[vname] = line.split()[col-1]
	return dictionary

def is_a_number(string):
	"""Takes a string and returns a boolean indicating whether it can be 
	converted to a float.
	"""
	try:
		float(string)
		return True
	except ValueError:
		return False

def merge_dicts(list_of_dicts):
    """Given any number of dicts, shallow copy and merge into a new dict,
	precedence goes to key value pairs in latter dicts.
	"""
    result = {}
    for dictionary in list_of_dicts:
        result.update(dictionary)
    return result

def print_ordered_dict(dictionary):
	"""Takes a dictionary and prints a sorted version of the dictionary.
	"""
	ordered_dictionary = collections.OrderedDict(sorted(dictionary.items()))
	for entry in ordered_dictionary:
		print(entry, ":", ordered_dictionary[entry])

def dict_intersect(dicts):
	"""Takes a list of dicts and returns a list with common labeled keys.
	"""
	common_keys = set([key for key in dicts[0].keys() if label in str(key)])
	for d in dicts[1:]:
		common_keys &= set([key for key in d.keys() if label in str(key)])
	result = sorted([key for key in common_keys])
	return result

def spartan_calculation(fpath):
	head, tail = os.path.split(fpath)  

	# if output file exists, calculation is not run
	if os.path.exists(fpath + ".spardir\\M0001\\output"): 
		print('Spartan calculation not run; output file already exists for ' + tail)
	else:
		# make file.spardir directory 
		if os.path.exists(fpath + ".spardir"):
			shutil.rmtree(fpath + ".spardir") 
		os.mkdir(fpath + ".spardir")

		# make empty _spartandir file 
		open(fpath + '.spardir\\_spartandir', 'w').close()

		# make empty M0001 directory
		os.mkdir(fpath + '.spardir\\M0001')

		# make empty _spartan file 
		open(fpath + '.spardir\\M0001\\_spartan', 'w').close()

		# extract molecular information from .spinput file created by user
		with open(fpath + '.spinput', 'r') as f:
			file_text = f.readlines() 
		molecular_data = in_between(file_text, 'M0001', 'ENDHESS')

		# make new input file in M0001 directory with calculation parameters
		with open(fpath + '.spardir\\M0001\\input', 'w') as f:
			f.writelines(header)
			f.writelines(molecular_data) 
			f.writelines(footer) 

		# run calculation via command line 
		if platform.system()=='Windows':
			# spartan_path = 'C:\\Program Files\\Wavefunction\\Spartan14v114'
			spartanx_path = spartan_path + '\\spartanx' 
			command = '"' + spartanx_path + '" --path "' + spartan_path + \
				'" --foreground-submit "' + fpath + '.spardir"'
			print('Calculation running for ' + tail + '...')
			print(command)
			call(command)
			print('Calculation complete')
		else:
			raise RuntimeError('Only Windows is currently supported')


def collinear(d, shared):
	"""Takes a dictionary containing atomic coordinates and a list of shared
	atoms and returns whether the labeled atoms are all collinear.
	"""
	if len(shared)<3:
		return True
	else:
		A = np.array(d[shared[0]])
		B = np.array(d[shared[1]])
		AB = np.linalg.norm(A-B)
		for i in range(2, len(shared)):
			C = np.array(d[shared[i]])
			AC = np.linalg.norm(A-C)
			BC = np.linalg.norm(B-C)
			# A, B, and C are collinear iff the largest of AB, AC, and BC 
			# is equal to the sum of the other two
			if AB+BC-AC<1e-8 or AB-BC+AC<1e-8 or -AB+BC+AC<1e-8:
				continue
			else:
				return False
		return True

def rigid_transform(A, B):
	"""Takes two matrices and returns rotation matrix and translation vector.
	
	Input: A, B
		A = Nx3 matrix of points
		B = Nx3 matrix of points
	Output: R, t
		R = 3x3 rotation matrix
		t = 3x1 translation vector

	>>> R, t = rigid_transform(A, B)
	>>> n = A.shape[0]
	>>> A_transform = np.dot(R, A.T) + np.tile(t, (1, n))
	>>> A_transform = A_transform.T

	A_transform should now be similar to B.  Note that since a geometry 
	optimization has been performed on the molecules, the coordinates will not 
	overlap exactly.
	"""
	assert A.shape==B.shape

	A = A.T
	B = B.T
	centroid_A = np.mean(A, axis=1, keepdims=True)
	centroid_B = np.mean(B, axis=1, keepdims=True)

	# center the points
	N = A.shape[1]
	AA = A - np.tile(centroid_A, (1, N))
	BB = B - np.tile(centroid_B, (1, N))

	# calculate covariance matrix
	H = np.dot(AA, BB.T)

	# perform singular value decomposition
	U, s, V = np.linalg.svd(H)
 	
	# check that SVD is correct
	S = np.zeros((3, 3))
	S[:3, :3] = np.diag(s)
	if not np.allclose(H, np.dot(U, np.dot(S, V))):
		raise ValueError('SVD failed')

	# calculate rotation matrix
	R = np.dot(V, U.T)

	# special reflection case
	if np.linalg.det(R) < 0:
		V[2,:] *= -1
		R = np.dot(V, U.T)

	# calculate translation vector
	t = np.dot(-R, centroid_A) + centroid_B

	return R, t

def compare_vibs(coords, vmodes, shared, i, j, varname, R2_min=0.50):
	"""Takes the indices of two molecules and extracts a correlation matrix of 
	the vibrational modes in those molecules.  Returns a list of the indices
	for vibrational modes that match.

	>>> result = compare_vibs(0, 1)
	>>> result
	[[5, 13], [6, 14], [7, 12], [9, 17], [10, 11], [12, 22], [13, 19]]

	In its current form, this function also saves a .csv file containing the 
	correlation matrix for heatmap generation in R.
	"""
	A_coords = np.array([coords[i][atom] for atom in shared])
	B_coords = np.array([coords[j][atom] for atom in shared])

	R, t = rigid_transform(A_coords, B_coords)

	R2 = []
	# for each combination of vibrational modes in the pair of molecules
	for m in range(len(vmodes[i])):
		for n in range(len(vmodes[j])):
			A_vibs = np.array([vmodes[i][m][atom] for atom in shared])
			B_vibs = np.array([vmodes[j][n][atom] for atom in shared])

			A_vibs_rot = np.dot(R, A_vibs.T)
			A_vibs_rot = A_vibs_rot.T

			A_vibs_rot_flat = [val for sublist in A_vibs_rot for val in sublist]
			B_vibs_flat = [val for sublist in B_vibs for val in sublist]

			# find pearson correlation for flattened coordinates
			pearson = pearsonr(A_vibs_rot_flat, B_vibs_flat)
			R2.append(pearson[0]**2)

	R2 = np.asarray(R2)
	R2 = np.reshape(R2, (len(vmodes[i]), len(vmodes[j])))
	# convert NaN's to 0
	R2[np.isnan(R2)] = 0

	# save R2 data to .csv
	plotname = 'R\\' + varname + '_corr' + str(j) + '.csv'
	np.savetxt(plotname, R2, delimiter=',')
	print(plotname, 'has been saved')

	# extract list of vibrational indices that match
	matched_vibs = []
	for row in range(R2.shape[0]):
		rowdata = R2[row, :]
		max_col = np.where(rowdata==max(rowdata))[0][0]
		coldata = R2[:, max_col]

		rowfreq = float(vmodes[i][row]['frequency'])
		colfreq = float(vmodes[j][max_col]['frequency'])

		# frequency calculations are not considered reliable below 500 cm-1
		if max(rowdata)==max(coldata) and max(rowdata)>=R2_min and \
		rowfreq>500 and colfreq>500:
			matched_vibs.append([row, max_col])

	row_labels = []
	col_labels = []
	for row in range(len(vmodes[i])):
		row_labels.append(float(vmodes[i][row]['frequency']))
	for col in range(len(vmodes[j])):
		col_labels.append(float(vmodes[j][col]['frequency']))

	plotname = 'R\\' + varname + '_corr' + str(j) + '_m1_labels.csv'
	with open(plotname, 'w') as output_file:
		writer = csv.writer(output_file, delimiter=',')
		writer.writerow(row_labels)

	plotname = 'R\\' + varname + '_corr' + str(j) + '_m2_labels.csv'
	with open(plotname, 'w') as output_file:
		writer = csv.writer(output_file, delimiter=',')
		writer.writerow(col_labels)
	return matched_vibs
		
def extract_coords(molecules):
	"""Takes a list of molecule names and returns a list of the atomic 
	coordinates and atomic weights for each molecule.
	"""
	coords = []
	for molecule in molecules:
		fpath = 'spartan_molecules\\' + molecule 

		# perform Spartan calculation
		spartan_calculation(fpath)
		
		with open(fpath + '.spardir\\M0001\\output', 'r') as f:
			output = f.readlines()

		# extract 3D atomic coordinates
		coords_text = in_between(output, 'Cartesian Coordinates', 'Point Group')
		coords_split = split(coords_text) # X,Y,Z coordinates of atoms

		# extract atomic coordinates
		acoords = {}
		for line in coords_split:
			if len(line)>0 and is_a_number(line[0]):
				acoords[line[2]] = np.array([float(coord) for coord in line[3:6]])
		coords.append(acoords)
		atoms = len(acoords)
	return coords



def count_atoms(molecule):
	"""Takes a molecule name and counts the number of atoms in the output file.
	"""
	fpath = 'spartan_molecules\\' + molecule
	with open(fpath + '.spardir\\M0001\\output', 'r') as f:
		output = f.readlines()

	# count number of atoms in molecule
	coords_text = in_between(output, 'Cartesian Coordinates', 'Point Group')
	coords_split = split(coords_text) 
	atoms = 0
	for line in coords_split:
		if len(line)>0 and is_a_number(line[0]):
			atoms += 1
	return atoms

def extract_atomic_weights(molecule):
	"""Takes a molecule name and creates a dictionary containing the atomic
	weight of each atom.
	"""
	fpath = 'spartan_molecules\\' + molecule
	with open(fpath + '.spardir\\M0001\\output', 'r') as f:
		output = f.readlines()

	# extract 3D atomic coordinates
	coords_text = in_between(output, 'Cartesian Coordinates', 'Point Group')
	coords_split = split(coords_text) 
	aw = {}
	for line in coords_split:
		if len(line)>0 and is_a_number(line[0]):
			aw[line[2]] = atomic_weights[line[1]]
	return aw

def get_shared_CH(molecules, shared):
	"""Takes a molecule name and creates a dictionary containing the types
	('H', 'C', etc.) of each atom.
	"""
	shared_CH = shared

	for molecule in molecules: 
		fpath = 'spartan_molecules\\' + molecule
		with open(fpath + '.spardir\\M0001\\output', 'r') as f:
			output = f.readlines()

		# extract atomic types
		coords_text = in_between(output, 'Cartesian Coordinates', 'Point Group')
		coords_split = split(coords_text) 

		labeled_CH = []
		for line in coords_split:
			if len(line)>0 and is_a_number(line[0]) and \
				line[2] in shared and (line[1]=='C' or line[1]=='H'):
					labeled_CH.append(line[2])
		shared_CH = list(set(shared_CH) & set(labeled_CH))
	return shared_CH


def extract_vmodes(molecules):
	"""Takes a list of molecule names and returns a list of dictionaries 
	containing the atomic movement, frequency, intensity, mode_id, and symmetry
	of each molecule.
	"""
	vmodes = []
	for molecule in molecules:
		fpath = 'spartan_molecules\\' + molecule

		with open(fpath + '.spardir\\M0001\\output', 'r') as f:
			output = f.readlines()

		# extract vibrational mode data (atomic coordinates)
		vibs_text = in_between(output, 'Normal Modes', 'Freq.') 
		vibs_split = split(vibs_text)

		# find [row, column] coords of 'X' where 'X','Y','Z' appears
		axes = ['X','Y','Z'] 
		x_coord = []
		row_num = 0
		for line in vibs_split: 
			row_num += 1
			for i in range(len(line)): 
				if line[i:i+len(axes)]==axes:
					x_coord.append([row_num, i])
		
		# extract table with vibrational mode data (mode, frequency, and intensity)
		# NOT the list of vibrational coordinates
		mode_data_text = in_between(output, 'Freq.', 'Standard')
		mode_data_split = split(mode_data_text)

		atoms = count_atoms(molecule)
		aw = extract_atomic_weights(molecule)

		# creates list of dictionaries containing vibrational modes
		# the xyz vector is weighted by the atomic weight of the atom
		mode_id = 0
		mol_vibs = []
		# for each vibrational mode
		for i in x_coord: 
			mode_id += 1
			row = i[0]
			col = i[1] + 1
			vib = {}
			for atom in range(atoms):
				xyz = []
				for dim in range(3): 
					aweight = aw[vibs_split[row+atom][0]]
					xyz.append(float(vibs_split[row+atom][col+dim]) * aweight)  
				vib[vibs_split[row+atom][0]] = xyz  # can do np.array() for clarity

			# add keys for symmetry, mode ID, frequency, and intensity
			vib['symmetry'] = vibs_split[row-2][int(col/3)]
			vib['mode_id'] = mode_id
			for line in mode_data_split:
				if len(line)>0 and is_a_number(line[0]) and int(line[0])==mode_id:
					vib['frequency'] = line[1]
					vib['intensity'] = line[2]
			mol_vibs.append(vib)
		vmodes.append(mol_vibs)
	return vmodes


def get_matched_vibs(coords, vmodes, shared, varname):
	"""Takes a list of molecules and returns a list of lists containing the 
	indices of matching vibrational modes.
	"""
	# compare first molecule with the rest successively, taking the intersection
	if len(vmodes)==1 or len(shared)<=1:
		matched_vibs = []
	else:
		matched_vibs = compare_vibs(coords, vmodes, shared, 0, 1, varname)
		if len(vmodes)>2:
			for i in range(2, len(vmodes)):
				new_vibs = compare_vibs(coords, vmodes, shared, 0, i, varname)
				for pair1 in matched_vibs:
					for pair2 in new_vibs:
						if pair1[0]==pair2[0]:
							pair1.append(pair2[1])
			matched_vibs = [vibs for vibs in matched_vibs if len(vibs)==len(vmodes)]
	return matched_vibs

def extract_descriptors(molecules, varname):
	"""Takes a list of molecule names and returns a list of dictionaries 
	containing molecular, atomic, and vibrational descriptors for each 
	molecule.
	"""
	coords = extract_coords(molecules)
	vmodes = extract_vmodes(molecules)
	shared = dict_intersect(coords)
	shared_CH = get_shared_CH(molecules, shared)
	matched_vibs = get_matched_vibs(coords, vmodes, shared, varname)

	descriptors = []
	mol_index = 0
	for molecule in molecules:
		fpath = 'spartan_molecules\\' + molecule
		with open(fpath + '.spardir\\M0001\\output', 'r') as f:
			output = f.readlines()

		# extract molecular descriptors
		keyword = 'Molecular volume:'
		vname = 'molecular_volume'
		molecular_volume = extract_line(output, keyword, vname, 3)

		keyword = 'Surface area:'
		vname = 'surface_area'
		surface_area = extract_line(output, keyword, vname, 3) 

		keyword = 'Ovality:'
		vname = 'ovality'
		ovality = extract_line(output, keyword, vname, 2)

		keyword = 'Atomic weight:'
		vname = 'molecular_weight'
		molecular_weight = extract_line(output, keyword, vname, 3)

		keyword = 'E(HOMO):'
		vname = 'E_HOMO'
		E_HOMO = extract_line(output, keyword, vname, 2)

		keyword = 'E(LUMO):'
		vname = 'E_LUMO'
		E_LUMO = extract_line(output, keyword, vname, 2)

		keyword = 'Electronegativity:'
		vname = 'electronegativity'
		electronegativity = extract_line(output, keyword, vname, 2)

		keyword = 'Hardness:'
		vname = 'hardness'
		hardness = extract_line(output, keyword, vname, 2) 

		# logP calculation failed for XPhos (and possibly other ligands)
		'''
		keyword = 'LogP (Ghose-Crippen):'
		vname = 'logP'
		logP = extract_line(output, keyword, vname, 3)
		'''

		keyword = 'Total Dipole:'
		vname = 'dipole_moment'
		dipole_moment = extract_line(output, keyword, vname, 3)

		# extract atomic descriptors for labeled atoms
		if len(shared)>0:
			charges_text = in_between(output, 'Atomic Charges:', 'Bond Orders')
			charges_split = split(charges_text)

			vname = 'electrostatic_charge'
			electrostatic_charge = extract_table(charges_split, shared, vname, 2, 4)

		else:
			electrostatic_charge = {}
			print('Atomic descriptors not extracted for', molecule)

		if len(shared_CH)>0:
			NMR_shift_text = in_between(output, 'NMR shifts', '<step')
			NMR_shift_split = split(NMR_shift_text)
			
			NMR_shift = extract_table(NMR_shift_split, shared_CH, 'NMR_shift', 2, 4)
		else:
			NMR_shift = {}
			print('NMR shifts not extracted for', molecule)

		# extract vibrational descriptors
		shared_vibs = {}
		if len(matched_vibs)==0:
			print('Vibrational descriptors not extracted for', molecule)
		else:
			for i in range(len(matched_vibs)):
				vib_index = matched_vibs[i][mol_index]
				label = 'V' + str(i+1)
				shared_vibs[label + '_frequency'] = \
					vmodes[mol_index][vib_index]['frequency']
				shared_vibs[label + '_intensity'] = \
					vmodes[mol_index][vib_index]['intensity']
		mol_index += 1

		# merge molecular, atomic, and vibrational descriptors
		mol_desc = merge_dicts([
			molecular_volume,
			surface_area,
			ovality,
			molecular_weight,
			E_HOMO,
			E_LUMO,
			electronegativity,
			hardness,
			# logP,
			dipole_moment,
			electrostatic_charge, 
			NMR_shift, 
			shared_vibs
		])
		descriptors.append(mol_desc)
	return descriptors

def var_dict(plates):
	"""Takes a list of Plate objects and returns a dictionary with variables as
	keys and lists of compounds as the values.  Gives a runtime error if the 
	variables are not consistent within or between plates.

	Example output:
	{'ArX': ['chlorobenzene', 'p-bromoacetophenone', 'p-chloroanisole'],
	'base': ['pyridine', 'triethylamine']}
	"""
	vars = {}
	v0 = set(plates[0].layout[0][0].conditions.keys())
	for plate in plates:
		for r in range(plate.rows):
			for c in range(plate.cols):
				# newkeys = set(plate.layout[r][c].conditions.keys())
				cond = plate.layout[r][c].conditions
				vn = set(cond.keys())
				if vn != v0:
					raise RuntimeError('Reaction variables not consistent ' \
						'between or within plates.')
				for key in cond.keys():
					if key not in vars.keys():
						vars[key] = [cond[key]]
					else:
						if cond[key] not in vars[key]:
							vars[key].append(cond[key])
	return vars

def get_descriptors(plates):
	"""Takes a list of Plate objects and returns a dictionary of the descriptors
	for every reagent in the plates.

	Example output:
	{'chlorobenzene': {descriptors}, 'pyridine': {descriptors}}
	"""
	vars = var_dict(plates)

	descriptors = {}
	# calculate descriptors for each variable
	for vname in vars.keys():
		d = extract_descriptors(vars[vname], vname)
		for i in range(len(vars[vname])):
			desc = d[i]
			newdesc = {}
			for key in desc.keys():
				newdesc[vname + '_' + key] = float(desc[key])
			descriptors[vars[vname][i]] = newdesc
	return descriptors

def export_reactions(plates):
	"""Takes a list of Plate objects and creates a .csv file with data where 
	each reaction is a row and the columns consist of the descriptors for every
	reagent in the reaction.
	"""
	descriptors = get_descriptors(plates)

	toCSV = []
	for plate in plates:
		for r in range(plate.rows):
			for c in range(plate.cols):
				cond = plate.layout[r][c].conditions
				rxn_desc = [descriptors[cond[key]] for key in cond.keys()]
				rxn_desc = merge_dicts(rxn_desc)
				toCSV.append(rxn_desc)


	keys = sorted(toCSV[0].keys())
	with open('R\\output_table.csv', 'w') as output_file:
		dict_writer = csv.DictWriter(output_file, keys)
		dict_writer.writeheader()
		dict_writer.writerows(toCSV)
	print('R\output_table.csv has been saved')

def export_for_pca(plates):
	"""Takes a list of Plate objects and creates a list of .csv files with the
	parameters for each variable type.
	"""
	vars = var_dict(plates)
	descriptors = get_descriptors(plates)

	# for variable in vars, create a table of descriptors
	for var in vars.keys():
		toCSV = []
		for mol in vars[var]:
			d = descriptors[mol]
			d['name'] = mol  # add molecule name to dictionary
			toCSV.append(d)

		keys = sorted(toCSV[0].keys())
		with open('R\\' + var + '.csv', 'w') as output_file:
			dict_writer = csv.DictWriter(output_file, keys)
			dict_writer.writeheader()
			dict_writer.writerows(toCSV)
		print('R\\' + var + '.csv has been saved')