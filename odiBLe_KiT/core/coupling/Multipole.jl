##################################################
# Structure used to compute the bare susceptibility
# coefficients
##################################################
struct structMultipole_type
    tabx ::Array{Float64,1} # Positions of the first particle
    tabxp::Array{Float64,1} # Positions of the second particle
    tabg ::Array{Float64,1} # Prefactors for the first particle
    tabgp::Array{Float64,1} # Prefactors for the second particle
    tabdeltaP0::Array{Float64,1} # Partial sum deltaP0
    tabdeltaP1::Array{Float64,1} # Partial sum deltaP1
    tabdeltaQ0::Array{Float64,1} # Partial sum deltaQ0
    tabdeltaQ1::Array{Float64,1} # Partial sum deltaQ1
    tabw ::Array{Int64,1} # Sorting of the arrays. ATTENTION, to the length. ATTENTION, indices are shifted by 1 because arrays starts at 1 in julia
end
##################################################
# Create a Multipole_type
##################################################
function structMultipole_create()
    return structMultipole_type(zeros(Float64,K_FT),zeros(Float64,K_FT), # tabx and tabxp
                                zeros(Float64,K_FT),zeros(Float64,K_FT), # tabg and tabgp
                                zeros(Float64,K_FT),zeros(Float64,K_FT), # tabdeltaP0 and tabdeltaP1
                                zeros(Float64,K_FT),zeros(Float64,K_FT), # tabdeltaQ0 and tabdeltaQ1
                                zeros(Int64,K_FT+2)) # tabw !! ATTENTION, the length is K_FT+2
end


##################################################
# Function to update the position and prefactor
# for a given particle
# @TODO
##################################################
function tabxtabg!(n::Int64,ra::Float64,
                   tabx::Array{Float64,1},
                   tabg::Array{Float64,1})
    #####
    Delta_u = 2.0/(K_FT) # Considered step distance
    Omega = Omegafromra(ra) # Value of the frequency

    ##### 
    # Warm-up for the first point
    #####
    u = -1.0 
    x = ra*fu(u) # Initial radius of the particle
    dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du 
    # Step 1
    k1 = 0.5*Delta_u*dtdu # Current velocity for theta
    # Step 2
    u += 0.25*Delta_u # Updating the time by half a timestep
    x = ra*fu(u) # Current location of the radius, r=r(u)
    dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du
    #####
    k2 = 0.5*Delta_u*dtdu # Current velocity for theta
    # Step 3
    k3 = k2 # Does not need to be updated
    ########################################
    # Step 4
    # First true point
    u += 0.25*Delta_u # Updating the time by half a timestep
    x = ra*fu(u) # Current location of the radius, r=r(u)
    dtdu = Omega * dthetadu_omega(ra,u) # Current value of dtheta1/du and dtheta2/du, always well-posed
    #####
    k4 = 0.5*Delta_u*dtdu # Current velocity for theta
    #####
    # Better guess for the value of theta(u=-1.0+Delta_u)
    # with a first RK4 warm-up integration for a small step
    theta = (k1 + 2.0*k2 + 2.0*k3 + k4)/(6.0) # Updating theta
    #####
    tabx[1] = x # First position
    tabg[1] = dtdu*cos(n*theta) # First prefactor
    #####
    for istep=2:K_FT # Loop over the steps of the integration !! ATTENTION, we only make (K_FT-1) steps
        #####
        # Integrating theta for one timestep
        ########################################
        # Step 1
        k1 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Step 2
        u += 0.5*Delta_u # Updating the time by half a timestep
        x = ra*fu(u) # Current location of the radius, r=r(u)
        dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du
        #####
        k2 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Step 3
        k3 = k2 # Does not need to be updated
        ########################################
        # Step 4
        u += 0.5*Delta_u # Updating the time by half a timestep
        x = ra*fu(u) # Current location of the radius, r=r(u)
        dtdu = Omega * dthetadu_omega(ra,u) # Current value of dtheta1/du and dtheta2/du, always well-posed
        #####
        k4 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Update
        theta += (k1 + 2.0*k2 + 2.0*k3 + k4)/(6.0) # Updating theta
        #####
        # Putting the values of the radius and prefactors
        tabx[istep] = x # Putting the current radius of the particle
        tabg[istep] = dtdu*cos(n*theta) # Putting the current prefactor of the particle
    end
end
##################################################
# Determines the order of the positions
# ATTENTION, this assumes that tabx and tabxp
# are already ordered
##################################################
function tabw!(tabw::Array{Int64,1},tabx::Array{Float64,1},tabxp::Array{Float64,1})
    tabw[1] = 0 # Initial term. !! ATTENTION, indices are shifted by one
    icount = 0 # Temporary counter for filling on tabw
    for j=1:K_FT # Filling in tabw. !! ATTENTION, tabR and tabRp have to be sorted
        while (icount < K_FT && tabx[icount+1] <= tabxp[j]) # !! ATTENTION, to the order of the tests to avoid segmentation faults
            icount += 1 # Updating the counter
        end
        tabw[j+1] = icount # ATTENTION, indices are shifted by one
    end
    tabw[K_FT+2] = K_FT # Last term. !! ATTENTION, indices are shifted by one
end
