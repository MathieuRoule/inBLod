include("./../Binning/bin.jl")
#####
# Bins list : center of the bin + width
# Different from ../Binning/hist_E.jl file 
#####
const NBINS_E = 25 									# Number of bins
const E_min = pot(0.0) 									# Minimal energy to consider
const E_max = pot(0.0) + 2.5 * G * MTOTBATH * LAMBDA 	# Maximal energy to consider
const dE_bin = (E_max-E_min)/NBINS_E 					# Bins width

const Ebin_list = zeros(Float64,NBINS_E) 				# Array of the bin centers
for k = 1:NBINS_E
    Ebin_list[k] = bincenter(k,E_min,dE_bin)
end

#####
# Arrays of particles energies
#####
if COUPLING=="Landau"
    const NPART4DEE = NTESTPART # We only consider the massless test particles in case of Landau simulation
else
    const NPART4DEE = NPART # With collective effects, all particle energies are considered
end

const tabE_binning = zeros(Float64,NPART4DEE) # Array of the energy of considered particles
const tabE_init_binning = zeros(Float64,NPART4DEE) # Array of the intial energy of considered particles

const whichbinE = zeros(Int64,NPART4DEE) # Bin associated to each particle
const bincount = zeros(Int64,NBINS_E) # Number of particle in each bin
const meanEinbin = zeros(Float64,NBINS_E) # Mean value of the energy in the bin

const meanDE2 = zeros(Float64,NBINS_E,NTESTSAVES) # Estimated mean on a bin on 1 realization