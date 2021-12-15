# ----------------------------------------------------------------------------
# This code supplements the following paper: 
#
# "Design of enzymatic cascade reactors through multi-objective 
#                  dynamic optimization"
#
# Authors: Leandros Paschalidis, Barbara Beer, Samuel Sutiono, Volker Sieber, 
# Jakob Burger
#
# The paper was submited to: Biochemical Engineering Journal.
# ----------------------------------------------------------------------------
'''This code must be used together with Tables 5-12 in the supplementary 
material of our paper. You can use this code to simulate and plot the 
various Pareto-optimal process schedules calculated in our paper.'''  

import numpy as np
from scipy.integrate import odeint
import math as mt
import scipy.optimize as opt
import pylab

#------------------------------------------------------------------ Parameters 

# Maximum reaction rates (U/mg)
Vmax = {}
Vmax['UDH'] =     221.333
Vmax['GlucD'] =     8.876
Vmax['KdgD'] =      5.109
Vmax['KgsalDH'] =  40.701
Vmax['NOX'] =      16.409
# Kinetic parameters (mM)
Km = {}
Km['UDH',1] =       0.0780
Km['UDH',7] =       0.5884
Km['GlucD',3] =     0.2945
Km['KdgD',4] =      0.4652
Km['KgsalDH',5] =   0.4812
Km['KgsalDH',7] =   0.1760
Km['NOX',7] =       0.1420
Km['NOX',8] =       0.0050
Km['NOX',9] =       0.0045
# Molecular weights (mg/mmol)
MW = {}
MW['UDH'] =     31210
MW['GlucD'] =   51010
MW['KdgD'] =    34790
MW['KgsalDH'] = 57700
MW['NOX'] =     51940
# First order decay constant for NOX (min^(-1))
kNOX = 0.030 
# Rate constant for the glucaro-1,4-lactone opening (min^(-1))
kII = 0.013 
# Oxygen solubilty (mM)
S9_Star = 1*0.2099/0.774

#----------------------------------------------------------- Control variables
'''You can find the values of the control variables in Tables 5-12 of the 
supplementary material of our paper. Select the values of the control variables 
for the plot you wish to reproduce and add them below. Hit run and the 
corresponding process schedule will be simulated and plotted.'''  

# Running time (min)
t_f = 479
# Initial D-glucuronate titer (mΜ)
S1_initial = 84
# Initial NAD+ titer (mΜ)
S7_initial = 3
# UDH titer (μΜ)
E_UDH_initial = 0.35
# GlucD titer (μΜ)
E_GlucD_initial = 0.47
# KdgD titer (μΜ)
E_KdgD_initial = 1.15
# KgsalDH titer (μΜ)
E_KgsalDH_initial = 1.45
# Initial NOX titer (μΜ)
E_NOX_initial = 2.69
# Supplementation magnitude 1 (mM/min)
A_1 = 29
# Supplementation magnitude 2 (mM/min)
A_2 = 25
# Supplementation magnitude 3 (mM/min)
A_3 = 29
# Supplementation time 1 (min)
t_1 = 188
# Supplementation time 2 (min)
t_2 = 285
# Supplementation time 3 (min)
t_3 = 376
# Oxygen mass transfer coefficient (min^(-1))
kLa = 1.2 


#-----------------------------------------------------------------------------
def Reaction_rates(X,t):
    v           = {}
    v['I']      = (Vmax['UDH']*X[1]*X[7])/((Km['UDH',7]*X[7])+(Km['UDH',1]*X[1])+(X[1]*X[7]))*MW['UDH']*X[10]/1000000
    v['II']     = kII*X[2]
    v['III']    = (Vmax['GlucD']*X[3])/(Km['GlucD',3]+X[3])*MW['GlucD']*X[11]/1000000 
    v['IV']     = (Vmax['KdgD']*X[4])/(Km['KdgD',4]+X[4])*MW['KdgD']*X[12]/1000000
    v['V']      = (Vmax['KgsalDH']*X[5]*X[7])/((Km['KgsalDH',7]*X[7])+(Km['KgsalDH',5]*X[5])+(X[5]*X[7]))*MW['KgsalDH']*X[13]/1000000 
    v['VI']     = (Vmax['NOX']*X[8]*X[9])/(Km['NOX',8]+X[8]*(1+X[7]/Km['NOX',7])*(Km['NOX',9]+X[9]))*MW['NOX']*X[14]/1000000
    return(v)

def NOX_deactivation_and_supplementation(X,t):
    r           = {}
    r['_d^NOX'] = kNOX*X[14]
    r['_s^NOX'] = (A_1*(mt.tanh(100*t/t_f-100*t_1/t_f)-mt.tanh(100*t/t_f-(100*t_1/t_f+4/t_f)))+A_2*(mt.tanh(100*t/t_f-100*t_2/t_f)-mt.tanh(100*t/t_f-(100*t_2/t_f+4/t_f)))+A_3*(mt.tanh(100*t/t_f-100*t_3/t_f)-mt.tanh(100*t/t_f-(100*t_3/t_f+4/t_f))))
    return(r)

def Oxygen_transfer(X,t):
    N           = {}
    N['_O2']    = kLa*(S9_Star-X[9])
    return(N)

#-------------------------------------------------------DIFFERENTIAL_EQUATIONS

def Material_balances(X,t):        
    v = Reaction_rates(X,t)
    N = Oxygen_transfer(X,t)
    r = NOX_deactivation_and_supplementation(X,t)
    dS1dt = -v['I']
    dS2dt = +v['I']-v['II']
    dS3dt = +v['II']-v['III']
    dS4dt = +v['III']-v['IV']
    dS5dt = +v['IV']-v['V']
    dS6dt = +v['V']
    dS7dt = +v['VI']-v['I']-v['V']
    dS8dt = -dS7dt
    dS9dt = -1/2*v['VI']+N['_O2']
    dEUDHdt      = 0
    dEGlucDdt    = 0
    dEKdgDdt     = 0
    dEKgsalDHdt  = 0
    dENOXdt      = -r['_d^NOX']+r['_s^NOX']
    dXdt = [0, dS1dt, dS2dt, dS3dt, dS4dt, dS5dt, dS6dt, dS7dt, dS8dt, dS9dt, dEUDHdt, dEGlucDdt, dEKdgDdt, dEKgsalDHdt, dENOXdt]
    return(dXdt)

#-----------------------------------------------------------------------SOLVER
t = np.linspace(0,t_f,100)
X0 = [0, S1_initial, 0, 0, 0, 0, 0, S7_initial, 0, S9_Star, E_UDH_initial , E_GlucD_initial , E_KdgD_initial , E_KgsalDH_initial, E_NOX_initial]
X = odeint(Material_balances,X0,t)

#----------------------------------------------------------------------PLOTTING
import matplotlib.pyplot as plt
plt.figure(figsize=(9,8))
plt.subplot(311)
plt.plot(t,X[:,1] ,'b', label='S$_1$')
plt.plot(t,X[:,2] ,'c', label='S$_2$')
plt.plot(t,X[:,3], 'm', label='S$_3$')
plt.plot(t,X[:,4] ,'y', label='S$_4$')
plt.plot(t,X[:,5] , 'r', label='S$_5$')
plt.plot(t,X[:,6], 'k', label='S$_6$') 
plt.plot(t,X[:,7] , 'g', label='S$_7$')
plt.plot(t,X[:,8], '0.75', label='S$_8$')
plt.xlabel('$\it{t}$ / (min)', fontsize=17)
plt.ylabel('$\it{S}$$_i$ / (mM)', fontsize=17)
plt.grid(False)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)

plt.subplot(312)
plt.plot(t,X[:,10], 'b', label='E$^{UDH}$')
plt.plot(t,X[:,11], 'r', label='E$^{GlucD}$')
plt.plot(t,X[:,12], 'c', label='E$^{KdgD}$')
plt.plot(t,X[:,13], 'g', label='E$^{KgsalDH}$')
plt.plot(t,X[:,14], 'm', label='E$^{NOX}$')
plt.xlabel('$\it{t}$ / (min)', fontsize=17)  
plt.ylabel('$\it{E}$$^{j}$ / (μM)', fontsize=17)
plt.grid(False)
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)

plt.subplot(313)
plt.plot(t,X[:,9], 'b', label='S$_{9}$')
plt.xlabel('$\it{t}$ / (min)', fontsize=17)
plt.ylabel('$\it{S}$$_{9}$ / (mM)', fontsize=17)
plt.grid(False)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=0.5)  
     
pylab.savefig('ProcessScheduleSimulation.PNG')
pylab.savefig('ProcessScheduleSimulation.pdf')
plt.show()