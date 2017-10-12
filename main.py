from classes.plate import Plate
import setup

# This function can be used to create a list with all entries
# Example: L(7,12) produces [7,8,9,10,11,12]
def L(first, last):
	return [i for i in range(first, last+1)]


# PLATE 1 ======================================================================

# control reactions and reactions with additive 7 are not included
plate1 = Plate(24,45)

# additives
plate1.fillRow([1,4,7,10], 'additive', '5-phenylisoxazole') # additive 2
plate1.fillRow([2,5,8,11], 'additive', 'ethyl-3-methylisoxazole-5-carboxylate') # additive 4
plate1.fillRow([3,6,9,12], 'additive', 'ethyl-5-methylisoxazole-3-carboxylate') # additive 6
plate1.fillRow([13,16,19,22], 'additive', '4-phenylisoxazole') # additive 1
plate1.fillRow([14,17,20,23], 'additive', '3-phenylisoxazole') # additive 3
plate1.fillRow([15,18,21,24], 'additive', '3-methylisoxazole') # additive 5

# ligands
plate1.fillRow([1,2,3,13,14,15], 'ligand', 'XPhos')
plate1.fillRow([4,5,6,16,17,18], 'ligand', 't-BuXPhos')
plate1.fillRow([7,8,9,19,20,21], 'ligand', 't-BuBrettPhos')
plate1.fillRow([10,11,12,22,23,24], 'ligand', 'AdBrettPhos')

# aryl halides
plate1.fillColumn([1,16,31], 'aryl_halide', '1-chloro-4-(trifluoromethyl)benzene') # ArX 1
plate1.fillColumn([2,17,32], 'aryl_halide', '1-bromo-4-(trifluoromethyl)benzene') # ArX 2
plate1.fillColumn([3,18,33], 'aryl_halide', '1-iodo-4-(trifluoromethyl)benzene') # ArX 3
plate1.fillColumn([4,19,34], 'aryl_halide', '1-chloro-4-methoxybenzene') # ArX 4
plate1.fillColumn([5,20,35], 'aryl_halide', '1-bromo-4-methoxybenzene') # ArX 5
plate1.fillColumn([6,21,36], 'aryl_halide', '1-iodo-4-methoxybenzene') # ArX 6
plate1.fillColumn([7,22,37], 'aryl_halide', '1-chloro-4-ethylbenzene') # ArX 7
plate1.fillColumn([8,23,38], 'aryl_halide', '1-bromo-4-ethylbenzene') # ArX 8
plate1.fillColumn([9,24,39], 'aryl_halide', '1-ethyl-4-iodobenzene') # ArX 9
plate1.fillColumn([10,25,40], 'aryl_halide', '2-chloropyridine') # ArX 10
plate1.fillColumn([11,26,41], 'aryl_halide', '2-bromopyridine') # ArX 11
plate1.fillColumn([12,27,42], 'aryl_halide', '2-iodopyridine') # ArX 12
plate1.fillColumn([13,28,43], 'aryl_halide', '3-chloropyridine') # ArX 13
plate1.fillColumn([14,29,44], 'aryl_halide', '3-bromopyridine') # ArX 14
plate1.fillColumn([15,30,45], 'aryl_halide', '3-iodopyridine') # ArX 15

# bases
plate1.fillColumn(L(1,15), 'base', 'P2Et')
plate1.fillColumn(L(16,30), 'base', 'BTMG')
plate1.fillColumn(L(31,45), 'base', 'MTBD')

# uncomment to view layout for plate1
# plate1.printLayout()


# PLATE 2 ======================================================================

# control reactions are not included
plate2 = Plate(32,45)

# additives
plate2.fillRow([1,5,9,13], 'additive', '5-methylisoxazole') # additive 8
plate2.fillRow([2,6,10,14], 'additive', 'benzo[c]isoxazole') # additive 10
plate2.fillRow([3,7,11,15], 'additive', '3,5-dimethylisoxazole') # additive 12
plate2.fillRow([4,8,12,16], 'additive', 'methyl-isoxazole-5-carboxylate') # additive 14
plate2.fillRow([17,21,25,29], 'additive', 'ethyl-isoxazole-3-carboxylate') # additive 9
plate2.fillRow([18,22,26,30], 'additive', 'ethyl-5-methylisoxazole-4-carboxylate') # additive 11
plate2.fillRow([19,23,27,31], 'additive', 'ethyl-isoxazole-4-carboxylate') # additive 13
plate2.fillRow([20,24,28,32], 'additive', 'benzo[d]isoxazole') # additive 15

# ligands
plate2.fillRow([1,2,3,4,17,18,19,20], 'ligand', 'XPhos')
plate2.fillRow([5,6,7,8,21,22,23,24], 'ligand', 't-BuXPhos')
plate2.fillRow([9,10,11,12,25,26,27,28], 'ligand', 't-BuBrettPhos')
plate2.fillRow([13,14,15,16,29,30,31,32], 'ligand', 'AdBrettPhos')

# aryl halides
plate2.fillColumn([1,16,31], 'aryl_halide', '1-chloro-4-(trifluoromethyl)benzene') # ArX 1
plate2.fillColumn([2,17,32], 'aryl_halide', '1-bromo-4-(trifluoromethyl)benzene') # ArX 2
plate2.fillColumn([3,18,33], 'aryl_halide', '1-iodo-4-(trifluoromethyl)benzene') # ArX 3
plate2.fillColumn([4,19,34], 'aryl_halide', '1-chloro-4-methoxybenzene') # ArX 4
plate2.fillColumn([5,20,35], 'aryl_halide', '1-bromo-4-methoxybenzene') # ArX 5
plate2.fillColumn([6,21,36], 'aryl_halide', '1-iodo-4-methoxybenzene') # ArX 6
plate2.fillColumn([7,22,37], 'aryl_halide', '1-chloro-4-ethylbenzene') # ArX 7
plate2.fillColumn([8,23,38], 'aryl_halide', '1-bromo-4-ethylbenzene') # ArX 8
plate2.fillColumn([9,24,39], 'aryl_halide', '1-ethyl-4-iodobenzene') # ArX 9
plate2.fillColumn([10,25,40], 'aryl_halide', '2-chloropyridine') # ArX 10
plate2.fillColumn([11,26,41], 'aryl_halide', '2-bromopyridine') # ArX 11
plate2.fillColumn([12,27,42], 'aryl_halide', '2-iodopyridine') # ArX 12
plate2.fillColumn([13,28,43], 'aryl_halide', '3-chloropyridine') # ArX 13
plate2.fillColumn([14,29,44], 'aryl_halide', '3-bromopyridine') # ArX 14
plate2.fillColumn([15,30,45], 'aryl_halide', '3-iodopyridine') # ArX 15

# bases
plate2.fillColumn(L(1,15), 'base', 'P2Et')
plate2.fillColumn(L(16,30), 'base', 'BTMG')
plate2.fillColumn(L(31,45), 'base', 'MTBD')


# PLATE 3 ======================================================================

# control reactions are not included
plate3 = Plate(32,45)

# additives
plate3.fillRow([1,5,9,13], 'additive', 'ethyl-3-methoxyisoxazole-5-carboxylate') # additive 23
plate3.fillRow([2,6,10,14], 'additive', '3-methyl-5-phenylisoxazole') # additive 17
plate3.fillRow([3,7,11,15], 'additive', 'N,N-dibenzylisoxazol-3-amine') # additive 19
plate3.fillRow([4,8,12,16], 'additive', 'methyl-5-(furan-2-yl)isoxazole-3-carboxylate') # additive 21
plate3.fillRow([17,21,25,29], 'additive', '5-(2,6-difluorophenyl)isoxazole') # additive 16
plate3.fillRow([18,22,26,30], 'additive', 'N,N-dibenzylisoxazol-5-amine') # additive 18
plate3.fillRow([19,23,27,31], 'additive', '5-methyl-3-(1H-pyrrol-1-yl)isoxazole') # additive 20
plate3.fillRow([20,24,28,32], 'additive', 'methyl-5-(thiophen-2-yl)isoxazole-3-carboxylate') # additive 22

# ligands
plate3.fillRow([1,2,3,4,17,18,19,20], 'ligand', 'XPhos')
plate3.fillRow([5,6,7,8,21,22,23,24], 'ligand', 't-BuXPhos')
plate3.fillRow([9,10,11,12,25,26,27,28], 'ligand', 't-BuBrettPhos')
plate3.fillRow([13,14,15,16,29,30,31,32], 'ligand', 'AdBrettPhos')

# aryl halides
plate3.fillColumn([1,16,31], 'aryl_halide', '1-chloro-4-(trifluoromethyl)benzene') # ArX 1
plate3.fillColumn([2,17,32], 'aryl_halide', '1-bromo-4-(trifluoromethyl)benzene') # ArX 2
plate3.fillColumn([3,18,33], 'aryl_halide', '1-iodo-4-(trifluoromethyl)benzene') # ArX 3
plate3.fillColumn([4,19,34], 'aryl_halide', '1-chloro-4-methoxybenzene') # ArX 4
plate3.fillColumn([5,20,35], 'aryl_halide', '1-bromo-4-methoxybenzene') # ArX 5
plate3.fillColumn([6,21,36], 'aryl_halide', '1-iodo-4-methoxybenzene') # ArX 6
plate3.fillColumn([7,22,37], 'aryl_halide', '1-chloro-4-ethylbenzene') # ArX 7
plate3.fillColumn([8,23,38], 'aryl_halide', '1-bromo-4-ethylbenzene') # ArX 8
plate3.fillColumn([9,24,39], 'aryl_halide', '1-ethyl-4-iodobenzene') # ArX 9
plate3.fillColumn([10,25,40], 'aryl_halide', '2-chloropyridine') # ArX 10
plate3.fillColumn([11,26,41], 'aryl_halide', '2-bromopyridine') # ArX 11
plate3.fillColumn([12,27,42], 'aryl_halide', '2-iodopyridine') # ArX 12
plate3.fillColumn([13,28,43], 'aryl_halide', '3-chloropyridine') # ArX 13
plate3.fillColumn([14,29,44], 'aryl_halide', '3-bromopyridine') # ArX 14
plate3.fillColumn([15,30,45], 'aryl_halide', '3-iodopyridine') # ArX 15

# bases
plate3.fillColumn(L(1,15), 'base', 'P2Et')
plate3.fillColumn(L(16,30), 'base', 'BTMG')
plate3.fillColumn(L(31,45), 'base', 'MTBD')


# CREATE OUTPUT FILE ===========================================================

# creates an output.csv file in the rxnpredict folder
# the plates must be populated with the same rxn components (base, ligand, etc.)
setup.export_reactions([plate1,plate2,plate3])
setup.export_for_pca([plate1,plate2,plate3])