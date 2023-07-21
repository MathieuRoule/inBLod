#!/bin/bash

JULIA=julia

# Experiment choice (to comment/uncomment)
exp="Diffusion"                # Orbital diffusion illustrations

# Coupling choice (to comment/uncomment)
# COUPLING="Landau"     # Without collective effects
COUPLING="Balescu-Lenard"             # With collective effects


if test $exp == "Diffusion"
then
    CODE=./RunDiffusion.jl
    PLOTCODE=./plotDiffusion.jl
    cd ./Illustration
else
    echo "UNKNOWN EXPERIMENT"
    exit
fi

#####
# PARAMETERS OF THE SIMULATION
#####
NBATHPART=10000
NTESTPART=10

# scheme symplectic (default)
# order 2 (default)
DT=0.001

TWARM=0.0
TMAX=50.0
TSAVESCALE="linear"
TENSAVE=0.1
TTESTSAVE=0.1
TFULLSAVE=50.0

SEED=1 # Random seed

FBATHSAMPLING="GTE"
FTESTSAMPLING="circle"
MSAMPLING="testnomassp"

FILENAME="./../results/${exp}_${FBATHSAMPLING}_${COUPLING}"

#####
# Run
#####
# Running the computations
${JULIA} ${CODE} --nbathpart $NBATHPART --ntestpart $NTESTPART --dt $DT --twarm $TWARM --tmax $TMAX --tsavescale $TSAVESCALE --tensave $TENSAVE --ttestsave $TTESTSAVE --tfullsave $TFULLSAVE --seed $SEED --fbathsampling $FBATHSAMPLING --ftestsampling $FTESTSAMPLING --msampling $MSAMPLING --coupling $COUPLING  --filename $FILENAME
wait

# Plotting the results dumped in the HDF5 file
${JULIA} ${PLOTCODE} --filename $FILENAME
wait

if test $exp == "Diffusion"
then
    open "${FILENAME}_diffusion.pdf"
fi