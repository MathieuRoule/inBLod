# odiBLe_Sim READ-ME

odiBLe_Sim is a julia code dedicated to the N-body simulation of 1D self-gravitating systems taking (BL) or not (Landau) collective effects into account.

The "core" directory contains the main parts of the N-body simulation :
* The quasi-stationary states potential and sampling technique,
* The forces' computation,
* The leap-frog integrator
* A standard Run file with dumping function to specify according to the quantities of interest,
* A standard Sumup file to merge the results from different computational core (dumped in different file) in one file.
    
Each experiment from the "experiments" directory extract information from the N-body simulations :
* "Distribution" saves the positions and velocities of test particles and/or the full distribution at scheduled times,
* "Flux" saves energy histograms of the distribution at scheduled times. From theses times series of histograms, one can then extract the orbital flux (Figure 3 in the paper),
* "EnergyDiffusion" saves the mean orbital (energy) diffusion in given bins in energy at scheduled times. The obtained time series can then be averaged over realization and post-treated to extract the diffusion coefficient measures (Figure 2 in the paper),
* "Binning" can be used in different experiments to extract position, velocity, position-velocity or energy histograms across the simulation.
* "Illustration" illustrate the orbital diffusion of test particles.
    
The expected arguments for a simulation are given and explained in "core/Args.jl".

To get familiar with the code, try running
```
bash launch_exp.sh
```
at "experiments/" location.
You can then modify the parameters for the prediction directly in the 
bash file and run it again to test different couplings (with or 
without collective effects) for example.

This example only provides an "Illustration" experiment.
To obtain acceptable diffusion coefficients or flux measurement, it is necessary to run extensive simulation.

One can design its own experiments taking an existing one as example. The generic 
"core/Run.jl" file explains the functions to specify for an experiment. 

