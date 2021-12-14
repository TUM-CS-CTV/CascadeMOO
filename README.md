# CascadeMOO
A multi-objective optimization approach for the design of enzymatic cascade reactors. Simulation and optimization codes are provided. 

The codes provided here supplement the paper: Design of enzymatic cascade reactors through multi-objective dynamic optimization. 

Three pieces of code are provided. 
A can be used to simulate the results provided in Tables 5-12 of the supplementary information. 
To use code A download it, copy the values of the control variables you wish from our supplementary information you wish to simulate and run it in your selected python enviroment. 

Codes B and C can be used to find Pareto-optimal process schedules for the enzymatic reactor discussed in the paper. 
Save both codes in the same directory in your PC. 
Make sure that you have downloaded Pyomo, glpk and IPOPT 3.11.1. 
Run code B in your selected python enviroment (e.g. spyder) and a set of 7 Pareto-optimal solutions will be produced. 
You can recreate all results in the paper and more by varying the values of the parameters: Φ^EC, Φ^EC and k_La.
