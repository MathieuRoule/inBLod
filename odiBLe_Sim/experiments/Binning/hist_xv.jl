if !(@isdefined binjl)  # If it as never been included
    include("./bin.jl")   # Bins handling functions
end
###################################################
# Histogram in velocity
###################################################
#####
# Position/Velocity Bins list : center of the bin + width
# If not already defined (for exemple if hist_x.jl or 
# hist_v.jl is included before)
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
if !(@isdefined NBINS_v)
    const NBINS_v = 80 # Number of bins
    const v_min = -2.0 * SIGMA # Minimal velocity considered (in the DF reference frame)
    const v_max = 2.0 * SIGMA # Maximal velocity considered (in the DF reference frame)
    const dv_bin = (v_max-v_min)/NBINS_v # Width of a bin

    const vbin_list = zeros(Float64,NBINS_v) # Array of the bin centers
    for k = 1:NBINS_v
        vbin_list[k] = bincenter(k,v_min,dv_bin)
    end
end

@timeit to "allocation" const hist_xv = zeros(Int64,NBINS_x,NBINS_v,NTESTSAVES)     # (x,v) histogram along simulation

###################################################
# Adding each particle to the bins in which it is 
#
# Forcing_BL --> using the bath particles
# Forcing_Landau --> using the massless test particles
###################################################
#####
# Adding the particle p to the right bins
#####
function hist_xv_add_p!(n::Int64,p::Int64)
    xp = x_meanDF_RF(p) # Position in the mean DF reference frame 
    vp = v_meanDF_RF(p) # Velocity w.r.t. the mean DF reference frame
    #####
    # Bins associated to the particle
    ##### 
    whichbinx = whichbin(xp,x_min,dx_bin) # Could possibly be higher than NBINS_x or lower than 1
    whichbinv = whichbin(vp,v_min,dv_bin) # Could possibly be higher than NBINS_v or lower than 1
    #####
    # Adding to the right bins
    #####
    if ( (1 <= whichbinx <= NBINS_x) && (1 <= whichbinv <= NBINS_v) )
        hist_xv[whichbinx,whichbinv,n] += 1
    end
end
#####
# Loop over particles Balescu-Lenard case
#####
function hist_xv_calc_BL!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = (NTESTPART+1):NPART # For all bath particles
        hist_xv_add_p!(n,p)
    end
end
#####
# Same in Landau case
#####
function hist_xv_calc_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 1:NTESTPART # For all test particles
        hist_xv_add_p!(n,p)
    end
end
#####
# Same in Forcing-Landau case
# Here we assume that the first test particle does not have to 
# be accounted for.
#####
function hist_xv_calc_Forcing_Landau!(n::Int64)
    #####
    # Histogram calculation
    #####
    for p = 2:NTESTPART # For all massless test particles (only)
        hist_xv_add_p!(n,p)
    end
end

#####
# Landau or BL
#####
if ( (COUPLING == "Balescu-Lenard") || (COUPLING == "Forcing_BL") )
    const hist_xv_calc! = hist_xv_calc_BL!
elseif (COUPLING == "Landau")
    const hist_xv_calc! = hist_xv_calc_Landau!
elseif (COUPLING == "Forcing_Landau") 
	const hist_xv_calc! = hist_xv_calc_Forcing_Landau!
end

include("./IO_xv.jl")