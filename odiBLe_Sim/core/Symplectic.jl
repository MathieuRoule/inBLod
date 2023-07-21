# To perform symplectic integrations of the motion

##################################################
# Function to perform one second-order
# symplectic integration for a time DT
# @IMPROVE? -- If needed the loops could be made in parallel
##################################################
function integrate_DT_SECONDORDER!()
    ####################
    # Step 1
    ####################
    for n=1:NPART
        @inbounds tabx[n] += 0.5*DT*tabv[n] # Updating tabtheta
    end
    #####
    @timeit to "forces_comp" tabvdot!() # Computing the velocities
    #####
    for n=1:NPART
        @inbounds tabv[n] += DT*tabvdot[n] # Updating tabv
    end
    ####################
    # Step 2
    ####################
    for n=1:NPART
        @inbounds tabx[n] += 0.5*DT*tabv[n] # Updating tabtheta
    end
    @inbounds TIME[1] += DT
end
##################################################
# Integration for one timestep
##################################################
if ORDER == 2
    integrate_DT! = integrate_DT_SECONDORDER! # Second-order symplectic integrator
else
    println("ERROR: UNKNOWN ORDER") # Unknown integration order
    exit() # Quitting the program
end