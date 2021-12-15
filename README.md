# CascadeMOO
A multi-objective optimization approach for the design of enzymatic cascade reactors. 
Simulation and optimization codes are provided.

## Setup
Install [Anaconda](https://www.anaconda.com/products/individual) on your PC.
Open the Power shell prompt in anaconda navigator and run the following commands:
 
`conda install -c conda-forge pyomo`

`conda install -c conda-forge glpk`

`conda install -c conda-forge ipopt 3.11.1 `

## Run the simulation code
You can use this code to simulate the results of our paper. You can find the results 
in Tables 5-12 in the supplementary material of our paper. Save the simulation code i
n your pc. Open it with spyder. Substitute the values of the control variables for 
the process schedule you wish to simulate. 

## Run the optimization codes
Save both files in the same directory. To run the optimization codes open both
files in spyder and run the xx file. You can produce all optimization results in our
paper and more by varying the values of the parameters: Φ^EC, Φ^CC and k_La.

## Publications
When using this work, please cite our paper:

[1] Design of enzymatic cascade reactors through multi-objective dynamic optimization.
Leandros Paschalidis, Barbara Beer, Samuel Sutiono, Volker Sieber, Jakob Burger
Submited to the Biochemical Engineering Journal 
