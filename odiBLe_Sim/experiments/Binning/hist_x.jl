if !(@isdefined binjl)  # If it as never been included
    include("./bin.jl")   # Bins handling functions
end
###################################################
# Histogram in position
###################################################
#####
# Position Bins list : center of the bin + width
# If not already defined (for exemple if hist_xv.jl is 
# included before)
#####
if !(@isdefined NBINS_x)
    const NBINS_x = 100                         # Number of bins
    const x_min = -2.5 * LAMBDA                 # Minimal position considered (in the DF reference frame)
    const x_max = 2.5 * LAMBDA                  # Maximal position considered (in the DF reference frame)
    const dx_bin = (x_max-x_min)/NBINS_x        # Width of a bin

    const xbin_list = zeros(Float64,NBINS_x)    # Array of the bin centers
    for k = 1:NBINS_x
        xbin_list[k] = bincenter(k,x_min,dx_bin)
    end
end

@timeit to "allocation" const hist_x  = zeros(Int64,NBINS_x,NTESTSAVES)    # x histogram along simulation

###################################################
# Adding each particle to the bins in which it is 
#
# Forcing_BL --> using the bath particles
# Forcing_Landau --> using the massless test particles
###################################################
#####
# Adding the particle p to the right bins
#####
function hist_x_add_p!(n::Int64,p::Int64)
    xp = x_meanDF_RF(p) # Position in the mean DF reference frame
    #####
    # Bins associated to the particle
    ##### 
    whichbinx = whichbin(xp,x_min,dx_bin) # Could possibly be higher than NBINS or lower than 1
    #####
    # Adding to the right bin
    #####
    if (1 <= whichbinx <= NBINS_x)
        hist_x[whichbinx,n] += 1
    end
end
#####
# Loop over particles Balescu-Lenard case
#####
function hist_x_calc_BL!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = (NTESTPART+1):NPART # For all bath particles
        hist_x_add_p!(n,p)
    end
end
#####
# Same in Landau case
#####
function hist_x_calc_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 1:NTESTPART # For all test particles
        hist_x_add_p!(n,p)
    end
end
#####
# Same in Forcing-Landau case
# Here we assume that the first test particle does not have to 
# be accounted for.
#####
function hist_x_calc_Forcing_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 2:NTESTPART # For all massless test particles (only)
        hist_x_add_p!(n,p)
    end
end

#####
# Landau or BL
#####
if ( (COUPLING == "Balescu-Lenard") || (COUPLING == "Forcing_BL") )
    const hist_x_calc! = hist_x_calc_BL!
elseif (COUPLING == "Landau")
    const hist_x_calc! = hist_x_calc_Landau!
elseif (COUPLING == "Forcing_Landau") 
	const hist_x_calc! = hist_x_calc_Forcing_Landau!
end

include("./IO_x.jl")