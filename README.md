# Instructions for Using `rxnpredict`

1.	Download and install the following programs:
    - Spartan ’14 V1.1.4
    - Python 3 (The anaconda distribution is recommended, as it has packages required for the software to run:  Download at https://www.anaconda.com/download/)
    - R (Download at https://cran.r-project.org/mirrors.html – choose any mirror link)
    - R Studio (Download at https://www.rstudio.com/products/rstudio/download/)
    - Sublime Text 3 (Download at https://www.sublimetext.com/3)
2.	Add Anaconda as a PATH variable so that Python will execute the scripts within Sublime Text 3.  
    - Navigate to “This PC” in File Explorer
    - Right click "This PC" → Properties → Advanced system settings → Environment Variables...
    - In User variables, click "Path" variable → Edit → New → type "path\to\Anaconda3" (no quotes – e.g., C:\Users\Derek\Anaconda3)
3.	Go to https://github.com/doylelab/rxnpredict.  On right, click "Clone or Download" and then "Download Zip".  This will download a local copy of the repository (folder) to your computer.
4.	All of the molecules whose properties will be used for modeling must first be drawn in Spartan:
	- Use the Spartan GUI to draw the molecules, saving them in the spartan_molecules folder (within the rxnpredict folder). 
	- Be sure to label any shared atoms within a substrate class (ligand, base, etc.) with a "\*".  You can do so by right clicking an atom, then click "Properties".  Change the label text at the bottom of the dialog box (e.g., \*C1).
	- Save the molecules in both .spardir format (for future editing) and .spinput format (this is what the program uses).
5.  Modify the python scripts:
	- In setup.py, change the value of spartan_path (line 16 only) to the path of the Spartan14v114.exe file.  Be sure to use \\\\ between folder names.
	- In main.py, describe the 2D layout of the plates you have run (line 10 onwards). Helpful syntax:
		- `plate_name = Plate(x,y)` where x is the number of rows and y is the number of columns.
		- `plate_name.fillRow([list], 'substrate_class', 'molecule_name')` where 'molecule_name' corresponds to the name of that molecule's .spinput file.  Replace `fillRow` with `fillColumn` to populate columns instead.  Note: If your plate design does not conform to one molecule per row/column, you can modify the Plate's dimensions accordingly.  In an extreme case, an Nx1 plate can specify the components of each reaction individually.
		- After the plates are filled (all "cells" of the plates must be populated with the same kinds of substrate_class), insert the following lines (where the plate names match the ones you have created):
		```sh
		setup.export_reactions([plate1,plate2,plate3])
		setup.export_for_pca([plate1,plate2,plate3])
		```
	- Run main.py (ctrl + B in Sublime Text) to create R\output_table.csv.  This table consists of one row per reaction and one column per descriptor per reaction component.  For example, a 2000-reaction screen with 20 base descriptors, 50 ligand descriptors, and 30 additive descriptors would generate a table with 2000 rows and 100 columns in R\output_table.csv.
6.  Train model in R:
	- Save yield data in a single column in a file named yields.csv (code reactions without data as NA).
	- Open and run analysis_template.R.

# How the Program Works
1.  To be added.