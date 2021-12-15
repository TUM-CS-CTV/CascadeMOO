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
import pyomo.environ as pe
from pyomo.dae import *
import Optimization2_CascadeMOO as FfPF
import pylab
import matplotlib.pyplot as plt

SpaceTimeYield=[]
CofactorConsumption=[]
EnzymeConsumption=[] 
FinalTime=[]
InitialS1Concentration=[]
InitialS7Concentration=[]
EUDH=[]
EGlucD=[]
EKdgD=[]
EKgsalDH=[]
ENOX=[]
t1=[]
t2=[]
t3=[]
A1=[]
A2=[]
A3=[]
TotalEnzymeConcentration=[] 
info1=[]
info2=[]
O1s=[]
O2s=[]
O3s=[]



i = 1
w1 = 1

# a for loop for the epsilon-EC constraint 
for ECi in [0, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150]:
    # initial values for the control variables 
    if ECi == 0: 
        tfi = 300
        S1i = 0
        S7i = 0
        EUDHi = 0
        EGlucDi = 0
        EKdgDi = 0
        EKgsalDHi = 0
        ENOXi = 0
        tau1i = 40
        tau2i = 60 
        tau3i = 80
        A1i = 0
        A2i = 0
        A3i = 0
    elif ECi == 0.025: 
        tfi = 300
        S1i = 0
        S7i = 0
        EUDHi = 0
        EGlucDi = 0
        EKdgDi = 0
        EKgsalDHi = 0
        ENOXi = 0
        tau1i = 40
        tau2i = 60 
        tau3i = 80
        A1i = 0
        A2i = 0
        A3i = 0
    else: 
        # the control variables are initialized from the solutions of the previous run
        # in order to ensure faster conversion
        tfi = model.tf
        S1i = model.S1[model.tau.first()]
        S7i = model.S7[model.tau.first()]
        EUDHi = model.EUDH[model.tau.first()]
        EGlucDi = model.EGlucD[model.tau.first()]
        EKdgDi = model.EKdgD[model.tau.first()]
        EKgsalDHi = model.EKgsalDH[model.tau.first()]
        ENOXi = model.ENOX[model.tau.first()]
        tau1i = model.tau1
        tau2i = model.tau2 
        tau3i = model.tau3
        A1i = model.A1
        A2i = model.A2
        A3i = model.A3
    model = FfPF.FunctionforPF(ECi,tfi,S1i,S7i,EUDHi,EGlucDi,EKdgDi,EKgsalDHi,ENOXi,tau1i,tau2i,tau3i,A1i,A2i,A3i)
    # solver selection
    solver=pe.SolverFactory('ipopt')
    # maximum iteration count limit 
    solver.options['max_iter'] = 100000
    # selection of the acceptible tolerance to be stricter than default 
    solver.options['acceptable_tol'] = 10**(-10)
    results = solver.solve(model, tee=True)
    
    # saving the results of each optimization run inside the for loop 
    SpaceTimeYield.append(pe.value(model.OBJ[1]))   
    EnzymeConsumption.append(pe.value(model.EC))
    CofactorConsumption.append(pe.value(model.CC))
    FinalTime.append(pe.value(model.tf))
    InitialS1Concentration.append(pe.value(model.S1[model.tau.first()]))
    InitialS7Concentration.append(pe.value(model.S7[model.tau.first()]))
    EUDH.append(pe.value(model.EUDH[model.tau.first()]))
    EGlucD.append(pe.value(model.EGlucD[model.tau.first()]))
    EKdgD.append(pe.value(model.EKdgD[model.tau.first()]))
    EKgsalDH.append(pe.value(model.EKgsalDH[model.tau.first()]))
    ENOX.append(pe.value(model.ENOX[model.tau.first()]))
    A1c = pe.value(model.A1)/pe.value(model.tf)
    A2c = pe.value(model.A2)/pe.value(model.tf)
    A3c = pe.value(model.A3)/pe.value(model.tf)
    t1c = pe.value(model.tau1)*pe.value(model.tf)/100
    t2c = pe.value(model.tau2)*pe.value(model.tf)/100
    t3c = pe.value(model.tau3)*pe.value(model.tf)/100
    A1.append(A1c)
    A2.append(A2c)
    A3.append(A3c)
    t1.append(t1c)
    t2.append(t2c)
    t3.append(t3c)

    TotalEnzymeConcentration.append(pe.value(model.SumEnzymes))
    info1.append(results.solver.termination_condition)
    info2.append(results.solver.status)
    O1s.append(pe.value(model.O1[1])) 
    O2s.append(pe.value(model.O2[1])) 
    O3s.append(pe.value(model.O3[1])) 


    # printing the results of each optimization on an individual text file 
    print('tau', file = open("ParetoOptimalPoint{}.txt".format(w1), "w+"))
    print(list(model.tau), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S1[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S2[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S3[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S4', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S4[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S5', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S5[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S6', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S6[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S7', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S7[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S8', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S8[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S9', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.S9[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vI[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vII[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vIII[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v4', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vIV[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v5', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vV[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('v7', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.vVI[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('rdNOX', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.rdNOX[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('rsNOX', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.rsNOX[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('NO2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.NO2[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('Ep1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.EUDH[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('Ep3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.EGlucD[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('Ep4', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.EKdgD[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('Ep5', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.EKgsalDH[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('Ep6', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print([pe.value(model.ENOX[jo]) for jo in model.tau], file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('tf', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.tf), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('tau1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(t1c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('tau2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(t2c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('tau3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(t3c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('A1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(A1c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('A2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(A2c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('A3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(A3c, file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('SumEnzymes', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.SumEnzymes), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('OBJF', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.OBJ[1]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a")) 
    print('Yield', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.Yield), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S1[0]', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.S1[model.tau.first()]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('S7[0]', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.S7[model.tau.first()]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('O1', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.O1[1]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('O2', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.O2[1]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print('O3', file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    print(pe.value(model.O3[1]), file = open("ParetoOptimalPoint{}.txt".format(w1), "a"))
    
    # definition of the real time variable for plotting 
    tfp = pe.value(model.tf)
    realtime = [element * tfp for element in list(model.tau)]

    # plotting the process schedules     
    plt.figure(figsize=(9,8))
    plt.subplot(311)
    plt.plot(realtime, [pe.value(model.S1[x]) for x in model.tau], label = 'S1', c = 'b', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S2[x]) for x in model.tau], label = 'S2', c = 'c', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S3[x]) for x in model.tau], label = 'S3', c = 'm', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S4[x]) for x in model.tau], label = 'S4', c = 'y', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S5[x]) for x in model.tau], label = 'S5', c = 'r', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S6[x]) for x in model.tau], label = 'S6', c = 'k', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S7[x]) for x in model.tau], label = 'S7', c = 'g', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.S8[x]) for x in model.tau], label = 'S8', c = '0.75',  ls = '-', lw = '2')
    plt.xlabel('$\it{t}$ / (min)', fontsize=17)
    plt.ylabel('$\it{S}$$_i$ / (mM)', fontsize=17)
    plt.grid(False)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)

    plt.subplot(312)
    plt.plot(realtime, [pe.value(model.EUDH[x]) for x in model.tau], label = 'UDH', c = 'b', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.EGlucD[x]) for x in model.tau], label = 'GlucD', c = 'r', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.EKdgD[x]) for x in model.tau], label = 'KdgD', c = 'c', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.EKgsalDH[x]) for x in model.tau], label = 'KgsalDH', c = 'g', ls = '-', lw = '2')
    plt.plot(realtime, [pe.value(model.ENOX[x]) for x in model.tau], label = 'NOX', c = 'm', ls = '-', lw = '2')
    plt.xlabel('$\it{t}$ / (min)', fontsize=17)  
    plt.ylabel('$\it{E}$$^{\mathrm{j}}$ / (μM)', fontsize=17)
    plt.grid(False)
    plt.yticks(fontsize=17)
    plt.xticks(fontsize=17)

    plt.subplot(313)
    plt.plot(realtime, [pe.value(model.S9[x]) for x in model.tau], label = 'S9', c = 'b', ls = '-', lw = '2')
    plt.xlabel('$\it{t}$ / (min)', fontsize=17)
    plt.ylabel('$\it{S}$$_{9}$ / (mM)', fontsize=17)
    plt.grid(False)
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.01, hspace=0.5)
    
    plt.tight_layout()

    pylab.savefig('ParetoOptimalPoint{}.PNG'.format(i))
    pylab.savefig('ParetoOptimalPoint{}.PDF'.format(i))
    plt.show() 
    
    i = i + 1
    w1 = w1 + 1

# printing an overview of the results of all runs 
print('SpaceTimeYield', file = open("Overview.txt", "w+"))
print(SpaceTimeYield, file = open("Overview.txt", "a"))
print('EnzymeConsumption', file = open("Overview.txt", "a"))
print(EnzymeConsumption, file = open("Overview.txt", "a"))
print('CofactorConsumption', file = open("Overview.txt", "a"))
print(CofactorConsumption, file = open("Overview.txt", "a"))
print('FinalTime', file = open("Overview.txt", "a"))
print(FinalTime, file = open("Overview.txt", "a"))
print('InitialS1Concentration', file = open("Overview.txt", "a"))
print(InitialS1Concentration, file = open("Overview.txt", "a"))
print('InitialS7Concentration', file = open("Overview.txt", "a"))
print(InitialS7Concentration, file = open("Overview.txt", "a"))
print('EUDH', file = open("Overview.txt", "a"))
print(EUDH, file = open("Overview.txt", "a"))
print('EGlucD', file = open("Overview.txt", "a"))
print(EGlucD, file = open("Overview.txt", "a"))
print('EKdgD', file = open("Overview.txt", "a"))
print(EKdgD, file = open("Overview.txt", "a"))
print('EKgsalDH', file = open("Overview.txt", "a"))
print(EKgsalDH, file = open("Overview.txt", "a"))
print('ENOX', file = open("Overview.txt", "a"))
print(ENOX, file = open("Overview.txt", "a"))
print('A1', file = open("Overview.txt", "a"))
print(A1, file = open("Overview.txt", "a"))
print('A2', file = open("Overview.txt", "a"))
print(A2, file = open("Overview.txt", "a"))
print('A3', file = open("Overview.txt", "a"))
print(A3, file = open("Overview.txt", "a"))
print('t1', file = open("Overview.txt", "a"))
print(t1, file = open("Overview.txt", "a"))
print('t2', file = open("Overview.txt", "a"))
print(t2, file = open("Overview.txt", "a"))
print('t3', file = open("Overview.txt", "a"))
print(t3, file = open("Overview.txt", "a"))

print('TotalEnzymeConcentration', file = open("Overview.txt", "a"))
print(TotalEnzymeConcentration, file = open("Overview.txt", "a"))
print('info1', file = open("Overview.txt", "a"))
print(info1, file = open("Overview.txt", "a"))
print('info2', file = open("Overview.txt", "a"))
print(info2, file = open("Overview.txt", "a"))
print('O1s', file = open("Overview.txt", "a"))
print(O1s, file = open("Overview.txt", "a"))
print('O2s', file = open("Overview.txt", "a"))
print(O2s, file = open("Overview.txt", "a"))
print('O3s', file = open("Overview.txt", "a"))
print(O3s, file = open("Overview.txt", "a"))