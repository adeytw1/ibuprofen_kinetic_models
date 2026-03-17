'''
Docstring for vary_catalyst_hydrogenation_rxtr

This Python simulation emulates a model of a hydrogenation reactor in an ibuprofen production process as dictated by literature.
The model utilizes a differential equation solver to determine species concentrations over time and finds an absolute maximum of concentration
and the batch time of that concentration by varying catalyst concentration.

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
#Volume = 0.5 # L 500 cm3 reactor volume (from paper)
p_H2 = 80.0 #bar (From paper) partial pressure of H2, constantly fed
w_cat = 3 #grams

# Initial conditions
C_IBAP_0 = 5.4 #mol/L (From paper)
T_0 = 373 #K (From paper) 373 K
#C_cat = w_cat / Volume # gcat/L 
alpha = 4426 # g/mol Pd (From paper), deactivation parameter

# ODEs

def odes(T, C):

    C_IBAP, C_H2O, C_IBPE, C_IBEB, C_Olig, Volume = C

    C_cat = w_cat / Volume # gcat/

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

    return [dC_IBAP_dt, dC_H2O_dt, dC_IBPE_dt, dC_IBEB_dt, dC_Olig_dt, 0]

# Solve the ODEs
t_span = (0,2000)  # Time span for the simulation
t_eval = np.linspace(t_span[0], t_span[1], 1400)  # Time points to evaluate

vol_input = 0.5 

C_0 = [C_IBAP_0, 0, 0, 0, 0, vol_input]  # Initial conditions, last term is volume???

sol = sp.integrate.solve_ivp(odes, t_span, C_0, t_eval=t_eval)

C_IBAP, C_H2O, C_IBPE, C_IBEB, C_Olig, Vol = sol.y

C_IBPE_w_vol = []
volumes = np.arange(0.247, 0.248, 0.0001)  # 0.247 L to 0.248 L, note optimized around 0.2474 L

max_times = []
max_concs = []

for V in volumes:
    C_0 = [C_IBAP_0, 0, 0, 0, 0, V]
    sol = sp.integrate.solve_ivp(odes, t_span, C_0, t_eval=t_eval)
    C_IBAP, C_H2O, C_IBPE, C_IBEB, C_Olig, Vol = sol.y

    C_IBPE_w_vol.append(C_IBPE)

    # max for this volume
    idx_max = np.argmax(C_IBPE)
    t_max = sol.t[idx_max]
    c_max = C_IBPE[idx_max]

    max_times.append(t_max)
    max_concs.append(c_max)

    print(f"Catalyst concentration = {(w_cat/V):.4f} L: max C_IBPE = {c_max:.4f} mol/L at t = {t_max:.1f} s")

# Overall max across all volumes
overall_idx = np.argmax(max_concs)
overall_V = volumes[overall_idx]
overall_t = max_times[overall_idx]
overall_c = max_concs[overall_idx]

print("\nOverall maximum:")
print(f"Max C_IBPE = {overall_c:.4f} mol/L at t = {overall_t:.3f} s for catalyst concentration of = {(w_cat/overall_V):.4f} gcat/L")

reactor_volume = (496.12*1000*overall_t)/(178.27*overall_c*3600)
print(f"Reactor volume is {reactor_volume:.4f} Liters")

# Plotting 
plt.figure(figsize=(10, 6))
for V, C_IBPE in zip(volumes, C_IBPE_w_vol):
    plt.plot(t_eval, C_IBPE, label=f'Catalyst conc. = {(w_cat/V):.1f} gcat/L')
plt.xlabel('Time (s)')
plt.ylabel('C_IBPE (mol/L)')
plt.title('IBPE Concentration vs Time for Different Reactor Volumes')
plt.legend()
plt.show()

'''
Best case scanrio! C_cat = 3/.2474 = 12.1261115602 gcat/L, produces 4.3169 mol/L IBPE in 1322.353 s, must make 496.12 kg IBPE per hour, therefore
need 234.82 L reactor.
'''


'''
Addressing the question
- Important differences from SuperPro. This model also models various side-reactions, the catalyst plays a role in the reaction,
what else about the reaction, what else?
- Batch reactors could be incorporated into a continuous-like system by arranging batch reactors in parallel such that at least one batch reactor
is expelling output while the other batch reactors are being loaded/processed. Essentially, if you had placed a black box over the system, it
would be indistinguishable from continuous.
'''
