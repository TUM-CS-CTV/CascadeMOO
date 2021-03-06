# CascadeMOO
A multi-objective optimization approach for the design of enzymatic cascade reactors. 
Simulation and optimization codes are provided.

## Setup
Install [Anaconda](https://www.anaconda.com/products/individual) in your PC.
PYOMO, GLPK and IPOPT 3.11.1 are required to run the optimization codes. These can be installed by running the 
following commands in the Powershell Prompt of anaconda navigator:
 
`conda install -c conda-forge pyomo`

`conda install -c conda-forge glpk`

`conda install -c conda-forge ipopt=3.11.1 `

## Run the simulation code
You can use the simulation code (Simulation_CascadeMOO.py) to reproduce all results in our paper. You can find 
the numerical values of the results in Tables 6-13 in the supplementary material of
our paper. Save the simulation code in your pc. Open it with Spyder and run it. You can vary the
values of the control variables (*t*<sub>f</sub>, *E*<sup>UDH</sup>, *E*<sup>GlucD</sup>, *E*<sup>KdgD</sup>, *E*<sup>KgsalDH</sup>, *E*<sup>NOX</sup>, *A*<sub>1</sub>, *A*<sub>2</sub>, *A*<sub>3</sub>, *t*<sub>1</sub>, *t*<sub>2</sub>, *t*<sub>3</sub>) and of the volumetric oxygen mass transfer coefficient (*k*<sub>L</sub>*a*) to simulate different process schedules. 

## Run the optimization codes
You can use the the optimization codes (Optimization1_CascadeMOO.py & Optimization2_CascadeMOO.py) to produce all optimization results in our
paper and more. Save both files in the same directory. Open both files in Spyder and 
run the Optimization1_CascadeMOO.py file. You can vary the values of the following parameters: *Φ*<sup>EC</sup>, *Φ*<sup>CC</sup> and 
*k*<sub>L</sub>*a* to produce different sets of Pareto-optimal solutions. 

## Publications
When using this work, please cite our paper:

Leandros Paschalidis, Barbara Beer, Samuel Sutiono, Volker Sieber, Jakob Burger,
Design of enzymatic cascade reactors through multi-objective dynamic optimization,
Biochemical Engineering Journal, 2022, 108384, ISSN 1369-703X, https://doi.org/10.1016/j.bej.2022.108384.
