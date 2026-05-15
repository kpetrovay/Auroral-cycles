This repository contains code and date to ensure reproducibility of the results
in the project
"Reconstruction of pretelescopic solar activity cycles from the auroral
record"
by Petrovay et al. (2026)

The code files are as follows:

runaursim.sh

    Bourne shell script to run aursim.py with a grid of input paramaters 
    F and M. The output is stored in mcsim_results.dat


aursim_methodtest.py

    Simulate a series of 1000 solar cycles and analyze the result by different
    procedures to identify cycle minima and other epochs.  Printed output
    provides hit rates to compare. Comment out the lines needed in the given
    procedure.

aursim.py

    Simulate a series of 1000 solar cycles and analyzes the result by the
    optimized procedure to identify cycle minima. 
    The command takes 2 parameters: F and M.
    The output is appended to the file mcsim_results.dat


mcsim_results.dat

    Tabulated output from runarsim.sh (via aursim.py). Columns are:
    F   M   1-year hit rate  2-year hit rate  F_emp  mean cycle length  st.dev.of cycle lengths
    

mccolmap_Femp.py

    Reads mcsim_results.dat and plots F_emp against F and M. 

mccolmap_percents_spline.py
	
    Reads mcsim_results.dat and plots the 2-year hit rate against F and M

aurora.py

    Analizes the auroral record in annual_aurorae_nokorean_1500_1650.xlsx by the
    optimized procedure to identify cycle minima. Output is printed and plotted.

aurrecon_data.zip

    Observational data used in the project. Upon unzipping into a folder,
    README_dat.txt describes the files. Copy or link the needed input file into
    the main folder containing the code. 

