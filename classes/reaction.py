class Reaction:

	def __init__(self):
		self.conditions = {}

	def addReagent(self, var, value):
		self.conditions[var] = value

	def printConditions(self):
		print(self.conditions)

	def getVariables(self):
		return self.conditions.keys()