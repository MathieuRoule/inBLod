# odiBLe_KiT READ-ME

odiBLe_KiT is a julia code dedicated to the computation of the inhomogeneous Balescu-Lenard equation of 1D self-gravitating systems (diffusion, friction, flux ...) taking (BL) or not (Landau) collective effects into account.

The "core" directory contains the main parts of the predictions computations :
    - The mean functions characterizing the quasi-stationary states,
    - Functions to compute the associated angle-action variables,
    - The bi-orthogonal basis and its Fourier transform,
    - The response matrix computation,
    - The coupling coefficients computations,
    - The diffusion, friction and flux coefficient computations.
    
Each experiment from the "experiments" directory extract different informations from the kinetic theory computations using the core code :
    - "BLpred" gives the total diffusion, friction and flux coefficients at given apocenters (orbits) as predicted, (Figures 2, 3 in the paper)
    - "Resonances_contribution" gives the contribution from each resonance (k,k') to the total diffusion, friction and flux at only one given apocenter. (Figures 6, 8 in the paper)
    
The expected arguments for a prediction are given and explained in "core/Args.jl".

To get familiar with the code, try running
```
bash launch_exp.sh
```
at "experiments/" location.

You can then modify the parameters for the prediction directly in the 
bash file and run it again to test different couplings (with or 
without collective effects) for example. You can also try other 
type of implemented experiments than diffusion prediction.
