from classes.reaction import Reaction
import setup
import csv

class Plate:

	def __init__(self, rows, cols):
		self.rows = rows
		self.cols = cols
		byrow = []
		for r in range(rows):
			bycol = []
			for c in range(cols):
				bycol.append(Reaction())
			byrow.append(bycol)
		self.layout = byrow

	def printLayout(self):
		print('Layout:', self.rows, 'x', self.cols)
		for row in range(len(self.layout)):
			for col in range(len(self.layout[0])):
				print('Row', row+1, 'Col', col+1, ':', end=' ')
				self.layout[row][col].printConditions()

	def fillBlock(self, row, col, var, value):
		for r in row:
			for c in col:
				self.layout[r-1][c-1].addReagent(var, value)

	def fillColumn(self, col, var, value):
		row = [i+1 for i in range(self.rows)]
		self.fillBlock(row, col, var, value)

	def fillRow(self, row, var, value):
		col = [i+1 for i in range(self.cols)]
		self.fillBlock(row, col, var, value)

	def fillTop(self, var, value):
		if self.rows % 2 == 0:
			row = [i+1 for i in range(int(self.rows/2))]
			col = [i+1 for i in range(self.cols)]
			self.fillBlock(row, col, var, value)
		else:
			raise RuntimeError('Cannot fill top half of plate because ' + \
				'there are an odd number of rows.')

	def fillBottom(self, var, value):
		if self.rows % 2 == 0:
			row = [i+1 for i in range(int(self.rows/2), self.rows)]
			col = [i+1 for i in range(self.cols)]
			self.fillBlock(row, col, var, value)
		else:
			raise RuntimeError('Cannot fill bottom half of plate because ' + \
				'there are an odd number of rows.')

	def fillLeft(self, var, value):
		if self.cols % 2 == 0:
			row = [i+1 for i in range(self.rows)]
			col = [i+1 for i in range(int(self.cols/2))]
			self.fillBlock(row, col, var, value)
		else:
			raise RuntimeError('Cannot fill left half of plate because ' + \
				'there are an odd number of columns.')
		
	def fillRight(self, var, value):
		if self.cols % 2 == 0:
			row = [i+1 for i in range(self.rows)]
			col = [i+1 for i in range(int(self.cols/2), self.cols)]
			self.fillBlock(row, col, var, value)
		else:
			raise RuntimeError('Cannot fill right half of plate because ' + \
				'there are an odd number of columns.')



	# OUT OF USE
	def sameVars(self):
		"""Returns whether a plate has the same variables in every reaction.
		"""
		vars = set(self.layout[0][0].conditions.keys())
		for row in self.layout:
			for rxn in row:
				if set(rxn.conditions.keys()) == vars:
					continue
				else:
					return False
		return True

	# OUT OF USE
	def getVarDict(self):
		"""Takes a plate and returns a dictionary containing a list of each
		reaction variable.
		
		Output:
		{'ArX': ['chlorobenzene', 'p-bromoacetophenone', 'p-chloroanisole'],
		'base': ['pyridine', 'triethylamine']}
		"""
		if self.sameVars()==False:
			raise RuntimeError('Not all reactions have the same variables.')

		vars = {}
		for r in range(self.rows):
			for c in range(self.cols):
				cond = self.layout[r][c].conditions
				for key in cond.keys():
					if key not in vars.keys():
						vars[key] = [cond[key]]
					else:
						if cond[key] not in vars[key]:
							vars[key].append(cond[key])
		return vars

	# OUT OF USE
	def getDescriptors(self):
		"""Takes a plate and returns a dictionary of the descriptors for every
		reagent in the plate.
		
		Output:
		{'chlorobenzene': {descriptors}, 'pyridine': {descriptors}}
		"""
		vars = self.getVarDict()

		descriptors = {}
		# calculate descriptors for each variable
		for vname in vars.keys():
			d = setup.extract_descriptors(vars[vname])

			for i in range(len(vars[vname])):
				desc = d[i]
				newdesc = {}
				for key in desc.keys():
					newdesc[vname + '_' + key] = float(desc[key])
				descriptors[vars[vname][i]] = newdesc

		# for key in descriptors.keys():
		# 	print('')
		# 	print(key)
		# 	setup.print_ordered_dict(descriptors[key])
		# 	print('')

		return descriptors

	# OUT OF USE
	def makeTable(self):
		"""Takes a Plate and creates a data table where each reaction is a 
		row and the columns consist of the descriptors for every reagent in 
		the reaction.
		"""
		descriptors = self.getDescriptors()

		toCSV = []
		for r in range(self.rows):
			for c in range(self.cols):
				cond = self.layout[r][c].conditions
				rxn_desc = [descriptors[cond[key]] for key in cond.keys()]
				rxn_desc = setup.merge_dicts(rxn_desc)
				toCSV.append(rxn_desc)

		keys = sorted(toCSV[0].keys())
		with open('R\\output_table.csv', 'w') as output_file:
			dict_writer = csv.DictWriter(output_file, keys)
			dict_writer.writeheader()
			dict_writer.writerows(toCSV)
		print('R\output_table.csv has been saved')