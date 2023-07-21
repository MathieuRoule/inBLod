##############################################################################################
# Force computation, different modes :
#   - Balescu-Lenard mode (bath and test particles both follows the instantaneous potential)
#   - Landau mode (bath particles follows the mean potential and test particles the instant.one)
#   - Forcing_BL (1st test particle follow the mean potential and all other particles the instant. one)
#   - Forcing_Landau (1st test particle and bath follow the mean potential and all other test particles the instant. one)
##############################################################################################

##################################################
# Instantaneous force 
# Caution, the different array have to be sorted/
# up-to-date at the time of the computation, specially
# the cumulative mass array
#
# Determining the force on all the particles
# Here, we denote as i, the index of the particle starting from the left
# The real index of that particle is then tabperm[i]=ind
# The INSTANTANEOUS force acting on that particle is
# vdot_ind = m*G*(N_right - N_left)
#          = m*G*(NPART - 2*N_left - 1)  [using N_right+N_left = N-1]
#          = m*G*(NPART - 2*i + 1)       [using N_left = i-1]
# With different masses, denoting m_ind the mass of the particle and 
# Mcum_i the cumulated mass which is equal to M_left + m_ind
# the force acting on this particle is vdot_ind*m_ind where
# vdot_ind = G*(M_right - M_left)
#          = G*(MTOT - 2*M_left - m_ind)     [using M_right+M_left = MTOT - m_ind]
#          = G*(MTOT - 2*Mcum_i + m_ind)     [using M_left = Mcum_i - m_ind]
##################################################
function instant_force(i::Int64)
    return @inbounds G*(MTOT[1] - 2.0*tabcumm[i] + tabm[tabperm[i]])
end

##################################################
# Mean Force
# Force from the mean field potential 
#
# Here, we denote as i, the index of the particle starting from the left
# The real index of that particle is then tabperm[i]=ind
# The INSTANTANEOUS force acting on that particle is
# vdot_ind = - dpotdx(x_ind)
# as x_ind is the position in the mean DF reference frame.
##################################################
function mean_force(ind::Int64)
    return - dpotdx(x_meanDF_RF(ind))
end

##################################################
# Sorting the particles from left to right and 
# Cumulated mass array actualization
##################################################
function sort_part_and_actualize_cumm!()
    # First, we find the appropriate sorting of the particles.
    # Here, tabperm[i] corresponds to the index of the particle
    # which is at the i^th position starting from the left
    sort!(tabperm,QuickSort,Perm(ForwardOrdering(),tabx)) # First, we sort the particles in-place.

    # The first particle on the left is treated separately for the cummulative mass array
    @inbounds ind_loc = tabperm[1] # Physical index of the particle that is that at the first position from the left
    @inbounds tabcumm[1] = tabm[ind_loc] # Cumulated mass for the first particle on the left

    for i=2:NPART
        @inbounds ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        @inbounds tabcumm[i] = tabcumm[i-1] + tabm[ind_loc] # Cumulated mass otherwise
    end
end




##################################################
# Balescu-Lenard mode :
# Function to compute the forces on all the particles
# This can be done efficiently simply by sorting
# all the particles.
# ATTENTION, this function should be checked again
# @IMPROVE? -- The filling-in of tabvdot is done
# in a non-sequential order: this likely can be done in a smarter way
# @IMPROVE? -- If needed, the loop could be made in parallel
##################################################
function tabvdot_BL!()
    # First we sort the particles and actualize the cumulative mass array
    sort_part_and_actualize_cumm!()

    for i=1:NPART # Loop over all the particles
        @inbounds ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left     
        @inbounds tabvdot[ind_loc] = instant_force(i) # All the particles follow the instantaneous potential
    end
end



############################################
# Landau mode :
# Same as Balescu-Lenard except that
# bath particles follow the mean potential
############################################

function tabvdot_Landau!()
    # First we sort the particles and actualize the cumulative mass array
    sort_part_and_actualize_cumm!()

    for i=1:NPART # Loop over the particles
        @inbounds ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        # If it is a test particle then it follows the instantaneous potential, else it follows the mean potential 
        # Caution : instantaneous force function needs the position from left i, where the mean one needs the 
        # index in the position/velocity arrays ind_loc=tabperm
        @inbounds tabvdot[ind_loc] = (ind_loc <= NTESTPART) ? instant_force(i) : mean_force(ind_loc)      
    end
end

############################################
# Forcing mode :
# The test particles follow the mean 
# potential and bath particles 
# follow the mean potential.
# Useful for a massive test particle 
# forcing an "external" perturbation on the 
# bath.
############################################

function tabvdot_Forcing_BL!()
    # First we sort the particles and actualize the cumulative mass array
    sort_part_and_actualize_cumm!()

    for i=1:NPART # Loop over the particles
        @inbounds ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        # If it is a test particle then it follows the mean potential, else it follows the instant. potential 
        @inbounds tabvdot[ind_loc] = (ind_loc <= NTESTPART) ? mean_force(ind_loc) : instant_force(i)
    end
end

############################################
# Forcing Landau mode :
# The first massive test particle and bath particles 
# follow the mean 
# potential and other test particles follow
# the instantaneous potential.  
# Useful for a massive test particle 
# forcing an "external" perturbation on the 
# bath without the bath seeing itself 
# (massless test particle see the inert bath + massive particle)
############################################

function tabvdot_Forcing_Landau!()
    # First we sort the particles and actualize the cumulative mass array
    sort_part_and_actualize_cumm!()

    for i=1:NPART # Loop over the particles
        @inbounds ind_loc = tabperm[i] # Physical index of the particle that is that at the i-th position from the left
        # If it is the first test particle or a bath particle then it follows the mean potential, 
        # else it follows the instant. potential 
        @inbounds tabvdot[ind_loc] = ((ind_loc == 1) || (ind_loc > NTESTPART)) ? mean_force(ind_loc) : instant_force(i)
    end
end

if COUPLING == "Balescu-Lenard"
    const tabvdot! = tabvdot_BL!
elseif COUPLING == "Landau"
    const tabvdot! = tabvdot_Landau!
elseif COUPLING == "Forcing_BL"
    const tabvdot! = tabvdot_Forcing_BL!
elseif COUPLING == "Forcing_Landau"
    const tabvdot! = tabvdot_Forcing_Landau!
end

