#!/bin/bash

# Number of available processors
NPROCS=4
export JULIA_NUM_THREADS=$NPROCS
export JULIA_CPU_THREADS=$NPROCS

JULIA=julia

# Experiment choice (to comment/uncomment)
exp="BLpred"                # Diffusion, friction, flux predictions
# exp="Rescontrib"          # Contribution of different resonances to the diffusion for a given orbit

if test $exp == "BLpred"
then
    CODE=./RunBL.jl
    cd ./BLpred
elif test $exp == "Rescontrib"
then
    CODE=./Rescontrib.jl
    cd ./Resonances_contribution
else
    echo "UNKNOWN EXPERIMENT"
    exit
fi


# Arguments for the prediction (some are useless for Landau prediction)
NPART=10000
NMAX=10

STATE="GTE"

MAPPING="polynomial"
KMF=100

# Coupling choice (to comment/uncomment)
COUPLING="Landau-Multipole"     # Without collective effects
# COUPLING="BL-PdS"             # With collective effects

KFT=1000

BASIS="CosSin"          # Useless for Landau-Multipole
NBASISELT=128           # Useless for Landau-Multipole
Lperiod=10.0            # Useless for Landau-Multipole

BLmethod="Legendre"     # Useless for Landau
NMAXRM=5                # Useless for Landau
ALEPH="neutral"         # Useless for Landau
BLramax=10.0            # Useless for Landau
KBL=50                  # Useless for Landau

FILENAME="./../results/${exp}_${STATE}_${COUPLING}"

# Running the computations
${JULIA} --threads $JULIA_NUM_THREADS ${CODE} --npart $NPART --nmax $NMAX --state $STATE --mapping $MAPPING --K_MF $KMF --coupling $COUPLING --K_FT $KFT --basis $BASIS --nbasiselt $NBASISELT --L_period $Lperiod --BLmethod $BLmethod --nmax_resmat $NMAXRM --aleph $ALEPH --BL_ra_max $BLramax --K_BL $KBL --filename $FILENAME
wait

# Plotting the results dumped in the HDF5 file
${JULIA} plot.jl --filename $FILENAME
wait


if test $exp == "BLpred"
then
    open "${FILENAME}_plot.pdf"
    open "${FILENAME}_plotDEE.pdf"
else
    open "${FILENAME}_plot.png"
fi

