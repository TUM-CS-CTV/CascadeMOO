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
from pyomo.environ import *
from pyomo.dae import *

def FunctionforPF(ECi,tfi,S1i,S7i,EUDHi,EGlucDi,EKdgDi,EKgsalDHi,ENOXi,
                  tau1i,tau2i,tau3i,A1i,A2i,A3i): 
    

    model = ConcreteModel()


    # total batch running time (min)
    model.tf = Var(bounds = (10,50000), initialize = tfi)     
    # scaled batch running time (min/min)                 
    model.tau = ContinuousSet(bounds=(0,1)) 
    # the reaction rates (mM/min)
    model.vI = Var(model.tau, within=NonNegativeReals)
    model.vII = Var(model.tau, within=NonNegativeReals)
    model.vIII = Var(model.tau, within=NonNegativeReals)
    model.vIV = Var(model.tau, within=NonNegativeReals)
    model.vV = Var(model.tau, within=NonNegativeReals)
    model.vVI = Var(model.tau, within=NonNegativeReals)
    # the NOX deactivation rate (mM/min)
    model.rdNOX = Var(model.tau, within=NonNegativeReals) 
    # the NOX supplementation rate (mM/min)
    model.rsNOX = Var(model.tau, within=NonNegativeReals)
    # the oxygen mass transfer rate (mM/min)
    model.NO2 = Var(model.tau, within=NonNegativeReals)
    # the substrate concentrations (mM)
    model.S1 = Var(model.tau, bounds = (0, 1000), initialize = S1i)
    model.S2 = Var(model.tau, bounds = (0, 1000))
    model.S3 = Var(model.tau, bounds = (0, 1000))
    model.S4 = Var(model.tau, bounds = (0, 1000))
    model.S5 = Var(model.tau, bounds = (0, 1000))
    model.S6 = Var(model.tau, bounds = (0, 1000))
    model.S7 = Var(model.tau, bounds = (0, 500), initialize = S7i)
    model.S8 = Var(model.tau, bounds = (0, 500))
    model.S9 = Var(model.tau, bounds = (0, 1000))
    # the enzyme concentrations (μM)
    model.EUDH= Var(model.tau, bounds = (0, 1000), initialize = EUDHi)
    model.EGlucD= Var(model.tau, bounds = (0, 1000), initialize = EGlucDi)
    model.EKdgD= Var(model.tau, bounds = (0, 1000), initialize = EKdgDi)
    model.EKgsalDH= Var(model.tau, bounds = (0, 1000), initialize = EKgsalDHi)
    model.ENOX= Var(model.tau, bounds = (0, 1000), initialize = ENOXi)
    # the first order derivatives of the substrate concentrations (mM/min)
    model.dS1dt = DerivativeVar(model.S1, wrt=model.tau, within=Reals)
    model.dS2dt = DerivativeVar(model.S2, wrt=model.tau, within=Reals)
    model.dS3dt = DerivativeVar(model.S3, wrt=model.tau, within=Reals)
    model.dS4dt = DerivativeVar(model.S4, wrt=model.tau, within=Reals)
    model.dS5dt = DerivativeVar(model.S5, wrt=model.tau, within=Reals)
    model.dS6dt = DerivativeVar(model.S6, wrt=model.tau, within=Reals)
    model.dS7dt = DerivativeVar(model.S7, wrt=model.tau, within=Reals)
    model.dS8dt = DerivativeVar(model.S8, wrt=model.tau, within=Reals)
    model.dS9dt = DerivativeVar(model.S9, wrt=model.tau, within=Reals)
    # the first order derivatives of the enzyme concentrations (μΜ)
    model.dEUDHdt = DerivativeVar(model.EUDH, wrt=model.tau, within=Reals)
    model.dEGlucDdt = DerivativeVar(model.EGlucD, wrt=model.tau, within=Reals)
    model.dEKdgDdt = DerivativeVar(model.EKdgD, wrt=model.tau, within=Reals)
    model.dEKgsalDHdt = DerivativeVar(model.EKgsalDH, wrt=model.tau, within=Reals)
    model.dENOXdt = DerivativeVar(model.ENOX, wrt=model.tau, within=Reals)
    # the supplementation rate for the first NOX supplementation (μΜ/min)
    model.NOXSup1 = Var(model.tau)
    # the supplementation rate for the second NOX supplementation (μΜ/min)
    model.NOXSup2 = Var(model.tau)
    # the supplementation rate for the third NOX supplementation (μΜ/min)
    model.NOXSup3 = Var(model.tau)
    # the total amount of NOX added  during the first supplementation (μΜ)
    model.O1 = Var(model.tau, within=NonNegativeReals)
    # the total amount of NOX added  during the second supplementation (μΜ)
    model.O2 = Var(model.tau, within=NonNegativeReals)
    # the total amount of NOX added  during the third supplementation (μΜ)
    model.O3 = Var(model.tau, within=NonNegativeReals)
    # the supplementation rate for the first NOX supplementation (μΜ/min)
    model.dO1 = DerivativeVar(model.O1, wrt=model.tau, within=Reals)
    # the supplementation rate for the second NOX supplementation (μΜ/min)
    model.dO2 = DerivativeVar(model.O2, wrt=model.tau, within=Reals)
    # the supplementation rate for the third NOX supplementation (μΜ/min)
    model.dO3 = DerivativeVar(model.O3, wrt=model.tau, within=Reals)
    # the scaled time of the first NOX supplementation (min/min)
    model.tau1 = Var(bounds = (0,100), initialize = tau1i)
    # the scaled time of the second NOX supplementation (min/min)
    model.tau2 = Var(bounds = (0,100), initialize = tau2i)
    # the scaled time of the third NOX supplementation (min/min)
    model.tau3 = Var(bounds = (0,100), initialize = tau3i)
    # the magnitude of the first NOX supplementation (μΜ/min)
    model.A1 = Var(within=NonNegativeReals, initialize = A1i)
    # the magnitude of the first NOX supplementation (μΜ/min)
    model.A2 = Var(within=NonNegativeReals, initialize = A2i)
    # the magnitude of the first NOX supplementation (μΜ/min)
    model.A3 = Var(within=NonNegativeReals, initialize = A3i)
    # the solubility of oxygen (mM)
    model.S9star = Var(within=NonNegativeReals)
    # the total concentration of all enzymes used during the batch (μΜ)
    model.SumEnzymes = Var(within=NonNegativeReals)
    # the yield (mM/mM)
    model.Yield = Var(within=NonNegativeReals)
    # the cofactor consumption (mM/min)
    model.CC = Var(within=NonNegativeReals)
    # the enzyme consumption (μΜ/min)
    model.EC = Var(within=NonNegativeReals, initialize = ECi)
    # the objective (Space-time yield) (mM/min)
    model.OBJ = Var(model.tau, within=NonNegativeReals)

    
    model.L = Set(initialize = ['UDH','GlucD','KdgD','KgsalDH','NOX'])
    model.M = Set(initialize = [1,2,3,4,5,6,7,8,9])

    # the molecular weights of all enzymes (mg/mM)
    mw = {}
    mw['UDH']    = 31210
    mw['GlucD']  = 51010
    mw['KdgD']   = 34790
    mw['KgsalDH'] = 57700
    mw['NOX']    = 51940
    model.mw = Param(model.L, initialize = mw)

    # the maximum reaction rates of all enzyme catalyzed reactions (U/mg)
    Vmax = {}
    Vmax['UDH'] = 221.331
    Vmax['GlucD'] = 8.876
    Vmax['KdgD'] = 5.109
    Vmax['KgsalDH'] = 40.701
    Vmax['NOX'] = 16.409
    model.Vmax = Param(model.L, initialize = Vmax)
    
    # the kinetic parameters of all enzyme catalyzed reactions (mM)
    Km = {}
    Km['UDH',1] = 0.0780
    Km['UDH',7] = 0.5884
    Km['GlucD',3] = 0.2945
    Km['KdgD',3] = 0.4652
    Km['KgsalDH',5] = 0.4812
    Km['KgsalDH',7] = 0.1760
    Km['NOX',7] = 0.1420
    Km['NOX',8] = 0.0050
    Km['NOX',9] = 0.0045
    model.Km = Param(model.L, model.M, initialize = Km)
    
    # the first order kinetic parameter for the lactone opening (min^(-1))
    model.kII = Param(initialize = 0.013)
    # the first order decay constant for the NOX deactivation (min^(-1))
    model.kNOX = Param(initialize = 0.03)
    # the epsilon constraint for the enzyme consumption (μΜ/min)
    model.ECc = Param(initialize = ECi)
    # the total pressure of the gas bubbles (atm)
    model.ptot = Param(initialize = 1)
    # the mole fraction of oxygen in the air bubbles (mol/mol)
    model.y9 = Param(initialize = 0.2099)
    # Henry's constant (m^3 atm / (mol))
    model.Hc = Param(initialize = 0.774) 
    # the volumetric mass transfer coefficient (min^(-1))
    model.kLa = Param(initialize = 1.2)
    
    # the reaction rate kinetics 
    def rr1(m, tau):
        return model.vI[tau] == model.EUDH[tau]*model.mw['UDH']*10**(-6)*(model.Vmax['UDH']*model.S1[tau]*model.S7[tau])/((model.Km['UDH',7]*model.S7[tau])+(model.Km['UDH',1]*model.S1[tau])+(model.S1[tau]*model.S7[tau]))
    model.rr1con = Constraint(model.tau, rule=rr1)
    
    def rr2(m, tau):
        return model.vII[tau] == model.kII*model.S2[tau]
    model.rr2con = Constraint(model.tau, rule=rr2)    
    
    def rr3(m, tau):
        return model.vIII[tau] == model.EGlucD[tau]*model.mw['GlucD']*10**(-6)*(model.Vmax['GlucD']*model.S3[tau])/(model.Km['GlucD',3]+model.S3[tau])
    model.rr3con = Constraint(model.tau, rule=rr3)
    
    def rr4(m, tau):   
        return model.vIV[tau] == model.EKdgD[tau]*model.mw['KdgD']*10**(-6)*(model.Vmax['KdgD']*model.S4[tau])/(model.Km['KdgD',3]+model.S4[tau])
    model.rr4con = Constraint(model.tau, rule=rr4)
    
    def rr5(m, tau): 
        return model.vV[tau] == model.EKgsalDH[tau]*model.mw['KgsalDH']*10**(-6)*(model.Vmax['KgsalDH']*model.S5[tau]*model.S7[tau])/((model.Km['KgsalDH',7]*model.S7[tau])+(model.Km['KgsalDH',5]*model.S5[tau])+(model.S5[tau]*model.S7[tau]))
    model.rr5con = Constraint(model.tau, rule=rr5)
    
    def rr6(m, tau):  
        return model.vVI[tau] == model.ENOX[tau]*model.mw['NOX']*10**(-6)*(model.Vmax['NOX']*model.S8[tau]*model.S9[tau])/((model.Km['NOX',8]+model.S8[tau]*(1+model.S7[tau]/model.Km['NOX',7]))*(model.Km['NOX',9]+model.S9[tau]))
    model.rr6con = Constraint(model.tau, rule=rr6)

    # the NOX deactivation and supplementation rates 
    def rE1(m, tau):
        return model.rdNOX[tau] == model.kNOX*model.ENOX[tau]
    model.rE1con = Constraint(model.tau, rule=rE1)
    
    def rE2(m, tau):
        return model.rsNOX[tau] == model.NOXSup1[tau]/model.tf + model.NOXSup2[tau]/model.tf + model.NOXSup3[tau]/model.tf
    model.rE2con = Constraint(model.tau, rule=rE2)
    
    def sup1NOX(m, tau):
        return model.NOXSup1[tau] == model.A1*((tanh(100*tau-model.tau1))-(tanh(100*tau-(model.tau1+4/model.tf))))
    model.sup1NOXcon = Constraint(model.tau, rule=sup1NOX)
    
    def sup1NOXo(m, tau):
        return model.dO1[tau] == model.NOXSup1[tau] 
    model.sup1NOXocon = Constraint(model.tau, rule=sup1NOXo)
    
    def sup2NOX(m, tau):
        return model.NOXSup2[tau] == model.A2*((tanh(100*tau-model.tau2))-(tanh(100*tau-(model.tau2+4/model.tf))))
    model.sup2NOXcon = Constraint(model.tau, rule=sup2NOX)
    
    def sup2NOXo(m, tau):
        return model.dO2[tau] == model.NOXSup2[tau] 
    model.sup2NOXocon = Constraint(model.tau, rule=sup2NOXo)
    
    def sup3NOX(m, tau):
        return model.NOXSup3[tau] == model.A3*((tanh(100*tau-model.tau3))-(tanh(100*tau-(model.tau3+4/model.tf))))
    model.sup3NOXcon = Constraint(model.tau, rule=sup3NOX)
    
    def sup3NOXo(m, tau):
        return model.dO3[tau] == model.NOXSup3[tau] 
    model.sup3NOXocon = Constraint(model.tau, rule=sup3NOXo)
        
    def ot1(m):
        return model.S9star == model.ptot*model.y9/model.Hc
    model.ot1con = Constraint(rule=ot1)
    
    def ot2(m, tau):
        return model.NO2[tau] == model.kLa*(model.S9star-model.S9[tau])
    model.ot2con = Constraint(model.tau, rule=ot2)    
    
    # the material balances for all substrates 
    def d1(m, tau):
        return model.dS1dt[tau] / model.tf == -model.vI[tau]
    model.d1con = Constraint(model.tau, rule=d1)
    
    def d2(m, tau):  
        return model.dS2dt[tau] / model.tf == +model.vI[tau]-model.vII[tau]
    model.d2con = Constraint(model.tau, rule=d2)
    
    def d3(m, tau):
        return model.dS3dt[tau] / model.tf == +model.vII[tau]-model.vIII[tau]
    model.d3con = Constraint(model.tau, rule=d3)
    
    def d4(m, tau):
        return model.dS4dt[tau] / model.tf == +model.vIII[tau]-model.vIV[tau]
    model.d4con = Constraint(model.tau, rule=d4)
    
    def d5(m, tau):
        return model.dS5dt[tau] / model.tf == +model.vIV[tau]-model.vV[tau]
    model.d5con = Constraint(model.tau, rule=d5)
    
    def d6(m, tau):
        return model.dS6dt[tau] / model.tf == +model.vV[tau]
    model.d6con = Constraint(model.tau, rule=d6)
    
    def d7(m, tau):
        return model.dS7dt[tau] / model.tf == +model.vVI[tau]-model.vI[tau]-model.vV[tau]
    model.d7con = Constraint(model.tau, rule=d7)
    
    def d8(m, tau):
        return model.dS8dt[tau] / model.tf == -(+model.vVI[tau]-model.vI[tau]-model.vV[tau])
    model.d8con = Constraint(model.tau, rule=d8)
    
    def d9(m, tau):
        return model.dS9dt[tau] / model.tf == +model.NO2[tau]-model.vVI[tau]/2
    model.d9con = Constraint(model.tau, rule=d9)
    
    # the material balances for all enzymes 
    def e1(m, tau):
        return model.dEUDHdt[tau]  == 0
    model.e1con = Constraint(model.tau, rule=e1)
    
    def e3(m, tau):
        return model.dEGlucDdt[tau] == 0
    model.e3con = Constraint(model.tau, rule=e3)
    
    def e4(m, tau):
        return model.dEKdgDdt[tau] == 0
    model.e4con = Constraint(model.tau, rule=e4)
    
    def e5(m, tau):
        return model.dEKgsalDHdt[tau] == 0
    model.e5con = Constraint(model.tau, rule=e5)

    def e6(m, tau):
        return model.dENOXdt[tau] / model.tf == -model.rdNOX[tau] +  model.rsNOX[tau]
    model.e6con = Constraint(model.tau, rule=e6)
    
    # calculation of the total enzyme concentration used 
    def c1(m):
        return model.SumEnzymes == model.EUDH[0] + model.EGlucD[0] + model.EKdgD[0] + model.EKgsalDH[0] + model.ENOX[0]  + model.O1[1] + model.O2[1] + model.O3[1] 
    model.c1con = Constraint(rule=c1)
        
    # definition of the yield 
    def c2(m):
        return model.Yield == (model.S6[1]-model.S6[0])/model.S1[0]
    model.c2con = Constraint(rule=c2)
    
    # a constraint on the yield 
    def c3(m):
        return model.Yield >= 0.95
    model.c3con = Constraint(rule=c3)
    
    # constraints to ensure that the supplementation times are ordered 
    def c4(m):
        return model.tau1 <= model.tau2
    model.c4con = Constraint(rule=c4)
    
    def c5(m):
        return model.tau2 <= model.tau3
    model.c5con = Constraint(rule=c5)
    
    # definition of the cofactor consumption 
    def c6(m):
        return model.CC == (model.S7[0] + model.S8[0]) / (model.tf+30)
    model.c6con = Constraint(rule=c6)
    
    # definition of the enzyme consumption
    def c7(m):
        return model.EC == model.SumEnzymes / (model.tf+30) 
    model.c7con = Constraint(rule=c7)
    
    # constraint on the cofactor consumption (changed manually)
    def c8(m):
        return model.CC <= 0.005
    model.c8con = Constraint(rule=c8)
    
    # constraint on the enzyme consumption (changed automatically)
    def c9(m, tau):
        return model.EC <= model.ECc
    model.c9con = Constraint(model.tau, rule=c9)
    
    # definition of the objective (space-time yield)
    def OBJ(m, tau):
        return model.OBJ[tau] == model.S6[1]/((model.tf+30))  
    model.OBJcon = Constraint(model.tau, rule=OBJ)
    
    # initial values for all substrates  
    model.ic = ConstraintList()
    model.ic.add(model.S2[model.tau.first()] == 0.001)
    model.ic.add(model.S3[model.tau.first()] == 0.001)
    model.ic.add(model.S4[model.tau.first()] == 0.001)
    model.ic.add(model.S5[model.tau.first()] == 0.001)
    model.ic.add(model.S6[model.tau.first()] == 0.001)
    model.ic.add(model.S8[model.tau.first()] == 0.001)
    model.ic.add(model.S9[model.tau.first()] == model.S9star)

    # selection of the objective for pyomo
    model.obj = Objective(expr=model.OBJ[1], sense=maximize)
    # selection of a discretization method
    discretizer = TransformationFactory('dae.finite_difference')
    # selection of the number of finite elements and discretization options 
    discretizer.apply_to(model, wrt=model.tau, nfe=100, scheme='BACKWARD')
   
    return model