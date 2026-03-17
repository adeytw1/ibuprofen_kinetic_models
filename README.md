## Ibuprofen Reactor Kinetic Models
This repository contains my kinetic models for three reactor present within our simulated ibuprofen plant. The three models are as follows.

# Semibatch Hydrogenation Reactor
_semibatch_time_hydrogenation_rxtr.py_
This is a semibatch hydrogenation reactor kinetic model. It is considered a semibatch reactor because hydrogen is constantly fed to maintain a constant hydrogen partial pressure. The hydrogenation step turns 4-isobutylacetophenone (IBAP) into the ibuprofen precursor 1-(4-isobutylphenyl)ethanol (IBPE). However, it also produces the unwanted byproducts oligomers (OLIG), 4-isobutylethylbenzene (IBEB) and water when IBAP reacts with catalyst in the presence of hydrogen gas through hydrogenolysis.

_vary_catalyst_hydrogenation_rxtr.py_
Additionally, a sensitivity analysis was performed regarding the concentration of catalyst. This was done to determine the best possible catalyst concentration numerically.

# Carbonylation Reactor
_carbonylation_rxtr.py_
This is a semibatch reactor kinetic model. It is considered a semibatch reactor because carbon monoxide is constantly fed to maintain a constant carbon monoxide partial pressure. The carbonylation step turns 1-(4-isobutylphenyl)ethanol (IBPE) into the desired product ibuprofen (IBN). IBPE reacts with free hydrogens from the acid to be converted into 5-isobutylstyrene (IBS), IBS is then convered to 1-(4-isobutylphenyl)ethyl chloride (IBPCl), which can freely convert back to IBS, IBPCl can then react wit hwater, carbon monoxide, and catalyst to product ibuprofen (IBN) and 3-(4-isobutylphenyl)propionic acid (IPPA).

For more information, please consult the final project deliverable.
