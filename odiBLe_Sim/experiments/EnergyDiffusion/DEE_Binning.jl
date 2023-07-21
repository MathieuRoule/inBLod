# Post-processing for DEE calculation

###################################################
# Calculation of individual particle energy
# in the mean field potential reference frame, then
# E_j = 1 / 2 * vp_j^2 + pot(xp_j).
# With vp_j = tabv[j] + vc and xp_j = tabx[j] + xc
# with (xc, vc) the offset of the simulation reference 
# frame to the mean field. (Already taken into account 
# in mean_energy(k::Int64) function)
# Input :
#   - tab -> to fill with the part. energy valeurs
###################################################

function part_energies!(tab::Array{Float64,1})
    for k = 1:NPART4DEE
        tab[k] = mean_energy(k)
    end
end

###################################################
# Initialization : initial energy and bin 
# associated to each particle
###################################################
if !(@isdefined binjl)  # If it as never been included
    include("./../Binning/bin.jl")   # Bins handling functions (whichbin function needed)
end
function init_binning!()
    #####
    part_energies!(tabE_init_binning) # Calculating the initial energy of each particle
    #####
    fill!(bincount,0.0) # Initialization of the bin count between different tries
    fill!(meanEinbin,0.0) # Initialization of the bin count between different tries
    for k = 1:NPART4DEE
        whichbinE[k] = whichbin(tabE_init_binning[k],E_min,dE_bin) # Could possibly be higher than NBINS_E --> to avoid in counting after
        if (1 <= whichbinE[k] <= NBINS_E)
            bincount[whichbinE[k]] += 1
            meanEinbin[whichbinE[k]] += tabE_init_binning[k]
        end
    end
    for p = 1:NBINS_E
        meanEinbin[p] /= (bincount[p] != 0) ? bincount[p] : 1.0 # Normalization of the mean
    end
end

###################################################
# Calculation of the mean value of
# (\Delta E)^2  = (E_j(t)-E_j(0))^2
# in each bin at a specific time given in 
# argument (index of save among the NPARTSAVES)
###################################################
function DE2_calc!(n::Int64)
    #####
    part_energies!(tabE_binning) # Energy calculation at current time
    #####
    # Mean estimation
    #####
    for k = 1:NBINS_E
        meanDE2[k,n] = 0.0 # Initialization
    end
    for k = 1:NPART4DEE
        bin = whichbinE[k]
        if (1 <= bin <= NBINS_E)
            meanDE2[bin,n] += (tabE_binning[k]-tabE_init_binning[k])^(2)
        end
    end
    for k = 1:NBINS_E
        meanDE2[k,n] /= (bincount[k] != 0) ? bincount[k] : 1.0 # Normalization of the mean
    end
end
