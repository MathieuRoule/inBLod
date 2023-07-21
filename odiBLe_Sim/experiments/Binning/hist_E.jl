if !(@isdefined binjl)  # If it as never been included
    include("./bin.jl")   # Bins handling functions
end
###################################################
# Histogram in energy
###################################################
#####
# Energy Bins list : center of the bin + width
#####
const NBINS_E = 30  									# Number of bins
const E_min = pot(0.0) 									# Minimal energy to consider (in the DF reference frame)
const E_max = pot(0.0) + 3.0 * G * MTOTBATH * LAMBDA 	# Maximal energy to consider (in the DF reference frame)
const dE_bin = (E_max-E_min)/NBINS_E 					# Bins width

const Ebin_list = zeros(Float64,NBINS_E) 				# Array of the bin centers
for k = 1:NBINS_E
    Ebin_list[k] = bincenter(k,E_min,dE_bin)
end

@timeit to "allocation" const hist_E  = zeros(Int64,NBINS_E,NTESTSAVES)     # Energy histogram along simulation

###################################################
# Adding each particle to the bins in which it is 
#
# Forcing_BL --> using the bath particles
# Forcing_Landau --> using the massless test particles
###################################################
#####
# Adding the particle p to the right bins
#####

function hist_E_add_p!(n::Int64,p::Int64)
    e_p = mean_energy(p) # Energy of the particle (according to the mean field)
    #####
    # Bins associated to the particle
    ##### 
    whichbine = whichbin(e_p,E_min,dE_bin) # Could possibly be higher than NBINS --> to avoid in counting after
    #####
    # Adding to the right bins
    #####
    if (1 <= whichbine <= NBINS_E)
        hist_E[whichbine,n] += 1
    end
end
#####
# Loop over particles Balescu-Lenard case
#####
function hist_E_calc_BL!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = (NTESTPART+1):NPART # For all bath particles
        hist_E_add_p!(n,p)
    end
end
#####
# Same in Landau case
#####
function hist_E_calc_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 1:NTESTPART # For all test particles
        hist_E_add_p!(n,p)
    end
end
#####
# Same in Forcing-Landau case
# Here we assume that the first test particle does not have to 
# be accounted for.
#####
function hist_E_calc_Forcing_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 2:NTESTPART # For all massless test particles (only)
        hist_E_add_p!(n,p)
    end
end

#####
# Landau or BL
#####
if ( (COUPLING == "Balescu-Lenard") || (COUPLING == "Forcing_BL") )
    const hist_E_calc! = hist_E_calc_BL!
elseif (COUPLING == "Landau")
    const hist_E_calc! = hist_E_calc_Landau!
elseif (COUPLING == "Forcing_Landau") 
	const hist_E_calc! = hist_E_calc_Forcing_Landau!
end

include("./IO_E.jl")