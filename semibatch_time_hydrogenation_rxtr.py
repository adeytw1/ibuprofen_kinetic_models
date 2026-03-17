'''
Docstring for semibatch_time_hydrogenation_rxtr

This Python simulation emulates a model of a hydrogenation reactor in an ibuprofen production process as dictated by literature.
The model utilizes a differential equation solver to determine species concentrations over time and finds an absolute maximum of concentration
IBPE.

All constants, model parameters, and models were provided through literature.

Modelling kinetics and deactivation for the selective hydrogenation of an aromatic ketone over Pd/SiO2
Nakul Thakar, Rob J. Berger, Freek Kapteijn, and Jacob A. Moulijn

Dynamic Plantwide Modeling, Uncertainty, and Sensitivity Analysis of a Pharmaceutical Upstream Synthesis: Ibuprofen Case Study
Frederico C. C. Montes, Krist Gernaey, and Gürkan Sin

CHEMICAL PROCESS DESIGN - CHEN 4520
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Constants
R = 8.314  # J/(mol K)
Volume = 0.5 # L 500 cm3 reactor volume (from paper)
p_H2 = 80.0 #bar (From paper) partial pressure of H2, constantly fed
w_cat = 3 #grams

# Initial conditions
C_IBAP_0 = 5.4 #mol/L (From paper)
T_0 = 373 #K (From paper)
C_cat = w_cat / Volume # gcat/L 
alpha = 4426 # g/mol Pd (From paper), deactivation parameter

def odes(T, C):

    C_IBAP, C_H2O, C_IBPE, C_IBEB, C_Olig = C

    # a is tracking degradation of catalyst in reference to oligomer due to spoling of catalytic sites ( that is more oligomer means less useful catalyst!!!)
    # alpha is 4426 g/mol Pd (from paper)

    k_1 = 1.14
    k_2 = 0.095
    k_3 = 0.024
    K_IBAP = 76.4 # L/mol
    K_H = 141 #atm-1
    K_H2O = 529 #mol-1

    a = np.exp(-alpha * (C_Olig/C_cat))
    
    # Reaction rates
    r_1 = (C_cat * a * k_1 * C_IBAP * p_H2)/((1 + K_IBAP*C_IBAP + np.sqrt(K_H*p_H2) + K_H2O*C_H2O)**2)
    r_2 = (C_cat * a * k_2 * C_IBPE * p_H2)/((1 + K_IBAP*C_IBAP + np.sqrt(K_H*p_H2) + K_H2O*C_H2O)**2)
    r_3 = (C_cat * a * k_3 * C_IBAP**2)/((1 + K_IBAP*C_IBAP + np.sqrt(K_H*p_H2) + K_H2O*C_H2O)**2)

    #Species balances
    dC_IBAP_dt = -r_1
    dC_IBPE_dt = r_1 - r_2
    dC_IBEB_dt = r_2
    dC_H2O_dt = r_2
    dC_Olig_dt = r_3

    return [dC_IBAP_dt, dC_H2O_dt, dC_IBPE_dt, dC_IBEB_dt, dC_Olig_dt]

# Solve the ODEs
t_span = (0,15000)  # Time span for the simulation
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # Time points to evaluate

C_0 = [C_IBAP_0, 0, 0, 0, 0]  # Initial concentrations

sol = sp.integrate.solve_ivp(odes, t_span, C_0, t_eval=t_eval)

C_IBAP, C_H2O, C_IBPE, C_IBEB, C_Olig = sol.y

# Plotting results


idx = np.argmax(C_IBPE)
time_max_IBPE = sol.t[idx]
conc_max_IBPE = C_IBPE[idx]
print(f'Max C_IBPE of {conc_max_IBPE} mol/L at time {time_max_IBPE} s')

plt.figure(figsize=(10, 6))
plt.plot(sol.t, C_IBAP, label='C_IBAP', color='blue')
plt.plot(sol.t, C_IBPE, label='C_IBPE', color='orange')
plt.plot(sol.t, C_IBEB, label='C_IBEB', color='green')
#plt.plot(sol.t, C_H2O, label='C_H2O', color='red')
plt.plot(sol.t, C_Olig, label='C_Olig', color='purple')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mol/L)')
plt.title('Concentration Profiles in Hydrogenation Reactor')
plt.legend()
plt.show()

# Vary concentration of catalyst, find that catalyst concentration has an effect on max IBPE concentration and time to reach it,
# Then find the volume to get the kg/hr.

# Need to ultimately produce 496.12 kg IBPE per hour as specified by design.