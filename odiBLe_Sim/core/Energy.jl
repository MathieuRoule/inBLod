# For energy and momentum calculation

##################################################
# Sorting the positions for energy computation
# and completing the cummulative mass and mass*position array
##################################################
function sort_en!()
    # Cf tabvdot!() for explanation
    sortperm!(tabperm,tabx,alg=QuickSort)

    # First particle on the left treated separately
    ind_loc = tabperm[1] # Physical index of the particle that is that at the first position from the left
    tabcummx[1] = tabm[ind_loc] * tabx[ind_loc] # Cumulated mass*position for the first particle on the left
    tabcumm[1] = tabm[ind_loc] # Cumulated mass for the first particle on the left

    for i=2:NPART # Loop over the other particles
        ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        tabcummx[i] = tabcummx[i-1] + tabm[ind_loc] * tabx[ind_loc] # Cumulated mass*position
        tabcumm[i] = tabcumm[i-1] +  tabm[ind_loc] # Cumulated mass
    end
end

##################################################
# Function to compute the total energy
# This can be done efficiently doing it when the 
# particles are sorted and taking it into account.
##################################################
function energy()
    #####
    # ATTENTION : we consider here that the particles are sorted 
    # and tabperm[i] gives the current location index of the i-th 
    # particle that is at the i-th position from the left.
    # Then we begin by sorting tabx.
    #####
    @timeit to "sort_en" sort_en!()
    en = 0.0
    for i=1:NPART # Loop over all particles from left to right
        ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        
        ##### 
        # By definition the total POTENTIAL energy is 
        # E_p = G *sum_{i<j} m_i * m_j * abs(x_i - x_j).
        # Using the fact that x is sorted ( i < j implies x_i < x_j ), 
        # we re-write E_p as 
        # E_p = G * sum_{j=1...N} m_j * P_j
        # where 
        # P_j = sum_{i=1...j} m_i * (x_j - x_i)
        #     = x_j * sum_{i=1...j} m_i - sum_{i=1...j} m_i * x_i
        #     = x_j * tabcumm[j] - tabcummx[j]
        # where tabcum and tabcummx are calculated in the sort_en!() 
        # Then we have 
        #
        # E_p = G * sum_{j=1...N} m_j * ( x_j * tabcumm[j] - tabcummx[j])
        #
        #
        # We then add the total kinetic energy which is naturally computed 
        # E_k = sum_{j=1...N} 0.5 * m_j * v_j^2
        #
        # In fact, tabm, tabx and tabv are not really sorted with particles from 
        # left to right but tabm[tabperm], tabx[tabperm], tabv[tabperm] gives this.
        # Noting ind(j) the physical index of the particle that is that at the i-th position from the left, we have
        #
        # E = E_p + E_k
        #   = sum_{j=1...N} m_{ind(j)} * ( G * (x_{ind(j)} * tabcumm[j] - tabcummx[j]) + 0.5 * vp_{ind(j)}^2 )
        # with vp = v - v_mf[1] the velocity in the mean DF reference frame.
        #####
        en += tabm[ind_loc] * ( G * ( tabx[ind_loc] * tabcumm[i] - tabcummx[i] ) + 0.5 * (v_meanDF_RF(ind_loc))^(2) ) 
    end
    return en
end

##################################################
# Function to compute the total momentum
# in the mean DF reference frame
##################################################
function momentum()
    p = 0.0
    for i=1:NPART 
        p += tabm[i] * v_meanDF_RF(i)
    end
    return p
end

##################################################
# Function to save energy and momentum
##################################################

function save_en!(k::Int64)
    if k > NENSAVES
        return nothing
    end
    tabE[k] = energy()
    tabP[k] = momentum()
    #####
    # Saving time at which the computation have been made
    #####
    save_current_time!(k,ensaved_times)
end

##################################################
# Function to compute the mean energy of the 
# particle ind.
# Using vp = v_ind - vmf correspond to 
# the velocity in the mean field reference frame
##################################################
function mean_energy(ind::Int64)
    return 0.5 * (v_meanDF_RF(ind))^(2) + pot(x_meanDF_RF(ind))
end
