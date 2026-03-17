'''
This Python simulation emulates a model of a carbonylation reactor in an ibuprofen production process as dictated by literature.
The model utilizes a differential equation solver to determine species concentrations over time.

All constants, model parameters, and models were provided through literature.

Kinetic Modeling of Carbonylation of 1-(4-Isobutylphenyl)ethanol Using a Homogeneous PdCl2(PPh3)2/TsOH/LiCl Catalyst System
Abdul M. Seayad, Jayasree Seayad, Patrick L. Mills, and Raghunath V. Chaudhari

Dynamic Plantwide Modeling, Uncertainty, and Sensitivity Analysis of a Pharmaceutical Upstream Synthesis: Ibuprofen Case Study
Frederico C. C. Montes, Krist Gernaey, and Gürkan Sin

CHEMICAL PROCESS DESIGN - CHEN 4520
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Constants
R = 8.314  # J/(mol K)
p_Reactor = 5.4 #MPa (From project) total pressure in reactor
p_CO = p_Reactor # Assume that the partial pressure of CO is equal to the total pressure of the reactor

# Initial Conditions
C_IBPE_0 = 1.13 # Initial concentration of IBPE in mol/L from project
C_HCl = 0.24 # Initial concentration of HCl in mol/L from project
C_H2O = 2.67 # Initial concentration of H2O in mol/L from project
Temp = 418.15 # Reactor Temperature in K from project, 388 K
C_catalyst = 1.121e-3 # Catalyst concentration in kmol/m3 from project, also 1.121e-3 mol/L

#Kinetic constants from paper 2?
k_i_0 = 6.07 #L/(mol-s)
k_ii_0 = 1.38e4 #L/(mol-s)
k_neg_ii_0 = 1.52e2 #L/(mol-s)
k_iii_0 = 1.64e5 #L2.43/(mol2.43-s)

Ea_1 = 21.77e3
Ea_2 = 43.62e3
Ea_neg_2 = 34.72e3
Ea_3 = 33.86e3

k_i = k_i_0 * np.exp(-Ea_1/(R*Temp)) # L/(mol-s)
k_ii = k_ii_0 * np.exp(-Ea_2/(R*Temp)) # L/(mol-s)
k_neg_ii = k_neg_ii_0 * np.exp(-Ea_neg_2/(R*Temp))
k_iii = k_iii_0 * np.exp(-Ea_3/(R*Temp)) #L/(mol-s)

#k_i = 0.0070 # L/(mol-s)
#k_ii = 0.0130 # L/(mol-s)
#k_neg_ii = 0.0038 # L/(mol-s)
#k_iii = 0.0629 #L/(mol-s)
K = 26.57
m = 0.43

def odes(T, C):
    C_IBPE, C_IBS, C_IBPCl, C_Prod = C

    wIPBE_wOverall = (C_IBPE*178.27)/(C_IBPE*178.27 + C_IBS*160.26 + C_IBPCl*196.72 + C_H2O*18 + C_HCl*36.46 + C_Prod*(2*206.28)) #Is this right?
    He = 1.77e-5*Temp - 4.4e-3 * wIPBE_wOverall

    #He = 5.798e-3 #This is an approximation, you'll have to do a calc later
    C_CO = p_CO * He

    r_A = (k_iii * C_IBPCl * C_H2O * C_CO * C_catalyst ** m)/(1+K*C_IBPCl) #Is this 0.43 or some other m value?
    r_i = k_i * C_IBPE * C_HCl #C_HCl is C_H+?
    r_ii = k_ii * C_IBS * C_HCl * C_HCl #C_HCl is C_H+ and C_Cl?
    r_neg_ii = k_neg_ii * C_IBPCl

    dC_IBPE_dt = -r_i
    dC_IBS_dt = r_i - r_ii + r_neg_ii
    dC_IBPCl_dt = r_ii - r_neg_ii - r_A
    dC_Prod_dt = r_A

    return [dC_IBPE_dt, dC_IBS_dt, dC_IBPCl_dt, dC_Prod_dt]

t_span = (0,45*60)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

#C_IBPE, C_IBS, C_IBPCl, C_Pl, C_H2O, C_HCl = C

'''
# Initial Conditions
C_IBPE_0 = 1.13 # Initial concentration of IBPE in mol/L from project
C_HCl_0 = 0.24 # Initial concentration of HCl in mol/L from project
C_H2O_0 = 2.67 # Initial concentration of H2O in mol/L from project
T = 388 # Reactor Temperature in K from project
C_catalyst = 1.121e-3 # Catalyst concentration in kmol/m3 from project, also 1.121e-3 mol/L
'''

C_0 = [C_IBPE_0, 0, 0, 0]

sol = sp.integrate.solve_ivp(odes, t_span, C_0, t_eval=t_eval)

C_IBPE, C_IBS, C_IBPCl, C_Prod = sol.y

reactor_size = (551.105*1000*max(t_span))/(206.29*C_Prod[-1]*3600)

# Prints summary of results
print("\nSUMMARY")
print(f"Time elapsed: {max(t_span)} seconds")
print(f"{C_Prod[-1]:.4f} mol/L of ibuprofen produced")
print(f"{C_IBPE[-1]:.4f} mol/L of IBPE produced")
print(f"{C_IBS[-1]:.4f} mol/L of IBS produced")
print(f"{C_IBPCl[-1]:.4f} mol/L of IBPCl produced")
print(f"Reactor size is {reactor_size:.4f} Liters")


# Plotting results
plt.figure(figsize=(10,6))
plt.plot(sol.t,C_IBPE,label="C_IBPE", color='blue')
plt.plot(sol.t,C_IBS,label="C_IBS", color='orange')
plt.plot(sol.t,C_IBPCl,label="C_IBPCl",color='red')
plt.plot(sol.t,C_Prod,label="C_Prod",color='purple')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mol/L)')
plt.title('Concentration Profiles')
plt.legend()
plt.show()

'''
Need to produce 551.105 kg/hr of ibuprofen by design specifications.
'''