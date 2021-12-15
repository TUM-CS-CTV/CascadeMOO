# CascadeMOO
A multi-objective optimization approach for the design of enzymatic cascade reactors. 
Simulation and optimization codes are provided.

Three pieces of code are provided. A can be used to simulate the results provided in 
Tables 5-12 of the supplementary information. B and C can be used to recreate our 
optimization results. 

# Setup
Install Anaconda on your PC
https://www.anaconda.com/products/individual
Open the Power shell prompt and run the following commands: 
conda install -c conda-forge pyomo
conda install -c conda-forge glpk
conda install -c conda-forge ipopt 3.11.1 

# Run the code
Save all files in the same directory. 
To run the simulation code, open it in spyder and run. 
You can simulate all results in the paper by selecting the optimal values for the control variables 
from Tables 5-12. To run the optimization codes open both files and run the xx file. You can recreate 
all optimization results by varying the values of the parameters: Φ^EC, Φ^EC and k_La.

# Publications
When using this work, please cite our publications:

[1] Design of enzymatic cascade reactors through multi-objective dynamic optimization.
Leandros Paschalidis, Barbara Beer, Samuel Sutiono, Volker Sieber, Jakob Burger
Biochemical Engineering Journal 2022.
