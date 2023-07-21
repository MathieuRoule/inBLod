##################################################
# Parameters for the integration
##################################################
@timeit to "allocation" const TIME = zeros(Float64,1) # Current time used in the integration
@timeit to "allocation" const MTOT = zeros(Float64,1) # Current total mass used in the integration

##################################################
# For recentering purpose
##################################################
const v_mf = zeros(Float64,1) # Momentum center of the mean DF in the simulation reference frame
# The velocity of a particle i in the mean DF reference is 
# vp = tabv[i] - v_mf[1]
# where tabv[i] is the velocity of this particle in the simulation reference frame.

##################################################
# Arrays describing the distribution of particles
##################################################
@timeit to "allocation" const tabm    = zeros(Float64,NPART) # Mass of the particle
@timeit to "allocation" const tabx    = zeros(Float64,NPART) # Position of the particles
@timeit to "allocation" const tabv    = zeros(Float64,NPART) # Velocities of the particles
@timeit to "allocation" const tabvdot = zeros(Float64,NPART) # Accelerations of the particles
@timeit to "allocation" const tabperm = zeros(Int64,NPART)   # Array used to store the particles' respective order
@timeit to "allocation" for i = 1:NPART
    tabperm[i] = i # Initialization of the permutation array so it is a permutation of 1:N from the beginning.
end

@timeit to "allocation" const tabcumm = zeros(Float64,NPART) # Cumulative mass of particles from left to right 
# the mass of the particle at i-th postition is included is the cumulated mass tabcummm[i]
# Then tabcumm[NPART] = MTOT[1]
@timeit to "allocation" const tabcummx = zeros(Float64,NPART) # Cumulative mass*position of particles from left to right 

##################################################
# For saving purpose
##################################################

#####
# Number and times of total energy and momentum saves at TENSAVE rate
# + arrays for total energy and momentum along simulation 
#####
if TSAVESCALE == "linear"
    const NENSAVES = floor(Int64,TMAX/TENSAVE)+1
elseif TSAVESCALE == "log2"
    const NENSAVES = max(0,floor(Int64,log2(TMAX/TENSAVE)))+2
elseif TSAVESCALE == "log10"
    const NENSAVES = max(0,floor(Int64,log10(TMAX/TENSAVE)))+2
end
@timeit to "allocation" const ensaved_times = zeros(Float64,NENSAVES)
@timeit to "allocation" const tabE    = zeros(Float64,NENSAVES) # Total energy along 1 simulation
@timeit to "allocation" const tabP    = zeros(Float64,NENSAVES) # Total momentum along 1 simulation

#####
# Number and times of test particles saves at TTESTSAVE rate
#####
if TSAVESCALE == "linear"
    const NTESTSAVES = floor(Int64,TMAX/TTESTSAVE)+1
elseif TSAVESCALE == "log2"
    const NTESTSAVES = max(0,floor(Int64,log2(TMAX/TTESTSAVE)))+2
elseif TSAVESCALE == "log10"
    const NTESTSAVES = max(0,floor(Int64,log10(TMAX/TTESTSAVE)))+2
end
@timeit to "allocation" const testsaved_times = zeros(Float64,NTESTSAVES)

#####
# Number and times of full particles saves at TFULLSAVE rate
#####
if TSAVESCALE == "linear"
    const NFULLSAVES = floor(Int64,TMAX/TFULLSAVE)+1
elseif TSAVESCALE == "log2"
    const NFULLSAVES = max(0,floor(Int64,log2(TMAX/TFULLSAVE)))+2
elseif TSAVESCALE == "log10"
    const NFULLSAVES = max(0,floor(Int64,log10(TMAX/TFULLSAVE)))+2
end

@timeit to "allocation" const fullsaved_times = zeros(Float64,NFULLSAVES)


