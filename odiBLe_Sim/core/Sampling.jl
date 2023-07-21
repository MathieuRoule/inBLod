##############################################################################################
# Centering the distribution / Change in reference frame
#
# To call after any addition of particle and before any integration.
#
# We do not change the position of any particle (so we do not move the mean-function distribution on the x-axis)
# but we center the discrete distribution in momentum to avoid a drift of the distribution.
# vc = sum(m_i*v_i) / sum(m_i) being the mean velocity of the discrete DF, we move every particle : v <- v - vc.
#
# Important assumption : 
# N is sufficiently large so the typical recentering in velocity is small enough 
# (comparing, e.g., the standard deviation dv of the recentering - which should be of mean almost 0 - 
# to the velocity bin size or dv^2/2 to the energy bin size).
# In this approx. the mean-field reference frame (in which we sample the particle) can be considered the same 
# as the simulation RF (tabx, tabv). 
# This assumption should be verified using bash_script/local/Sampling_test.sh script with the parameters 
# of the considered simulation. One can then verify that the sampling is done correctly through 
# the dumped histograms.
##############################################################################################
function centered_frame!()
    vc = 0.0 # Current total momentum
    for i=1:NPART 
        vc += tabm[i] * tabv[i]
    end
    vc /= MTOT[1] # Mean velocity (Total momentum / Total mass)

    for i=1:NPART # Move every particle
        tabv[i] -= vc # total momentum is now 0.
    end
    v_mf[1] -= vc # Actualizing the stored value of the mean DF momentum center in the simulation reference frame
end


#####################################
# Useful changes of coordinates
#####################################
#####
# Functions giving the position (resp. velocity) of the particle i in the mean DF reference frame
#####
function x_meanDF_RF(i::Int64)
    return tabx[i] 
end
function v_meanDF_RF(i::Int64)
    return tabv[i] # We here consider the RF to be the same
end
#####
# Function giving the position (resp. velocity) in the mean DF reference frame given the position 
# (resp. velocity) in the simulation reference frame.
# Therefore x/v_meanDF_RF(i) = x/v_from_sim_to_meanDF_RF(tabx/v[i]).
#####
function x_from_sim_to_meanDF_RF(x::Float64)
    return x
end
function v_from_sim_to_meanDF_RF(v::Float64)
    return v # We here consider the RF to be the same
end
#####
# Function giving the position (resp. velocity) in the simulation reference frame given the position 
# (resp. velocity) in the mean DF reference frame (usefull for sampling after a recentering).
#####
function x_from_meanDF_to_sim_RF(xp::Float64)
    return xp
end
function v_from_meanDF_to_sim_RF(vp::Float64)
    return vp # We here consider the RF to be the same
end

##############################################################################################
# Usefull functions (potential, derivative ....)
##############################################################################################

if (INITFBATH == "GTE") || (INITFTEST == "GTE")
    include("states/GTE.jl")
end
if (INITFBATH == "Plummer") || (INITFTEST == "Plummer")
    include("states/Plummer.jl")
end
if (INITFBATH == "Polytrope") || (INITFTEST == "Polytrope")
    include("states/Polytrope.jl")
end

##################################################

if INITFBATH == "GTE"
    const pot = potential_GTE
    const dpotdx = dpotdx_GTE
    const inv_pot = inv_pot_GTE
elseif INITFBATH == "Polytrope"
    const pot = potential_pol
    const dpotdx = dpotdx_pol
    const inv_pot = inv_pot_pol
elseif INITFBATH == "Plummer"
    const pot = potential_plum
    const dpotdx = dpotdx_plum
    const inv_pot = inv_pot_plum
elseif INITFBATH == "rand"
    const pot = potential_GTE # TO DO
    const dpotdx = dpotdx_GTE # TO DO
    const inv_pot = inv_pot_GTE # TO DO
elseif INITFBATH == "bar"
    const pot = potential_GTE # TO DO
    const dpotdx = dpotdx_GTE # TO DO
    const inv_pot = inv_pot_GTE # TO DO
end

##############################################################################################
# Initialization of the position and velocities
##############################################################################################

##################################################
# For bath particles
##################################################

#####
# Sampling uniformly in a square of the phase space
# such that total energy is 3/4 * G * (MTOTBATH)^(2) * LAMBDA
# as for GTE (see Violent_relaxation MMA Notebook)
#####
#####
# Sampling the particle p
#####
function sampling_rand_p!(p::Int64)
    tabx[p] = 1.75 * LAMBDA * (2.0 * rand() - 1.0)
    tabv[p] = v_from_meanDF_to_sim_RF(1.0 * SIGMA * (2.0 * rand() - 1.0))
end
#####
# Sampling the bath particles
#####
function sampling_rand_bath!()
    for p=(NTESTPART+1):NPART
        sampling_rand_p!(p)
    end
end

#####
# Sampling bar for phase mixing similar to wendy
# x in [-pi/2 , pi/2] uniform
# v = -V0 sin(x)
#####
#####
# Sampling the particle p
#####
function sampling_bar_p!(p::Int64,N::Int64)
    v0 = 0.001
    tabx[p] = ( (p-1)/(N-1) - 0.5) 
    tabv[p] = v_from_meanDF_to_sim_RF(-v0 * sin(tabx[p]))

end
#####
# Sampling the bath particles
#####
function sampling_bar_bath!()
    for p=(NTESTPART+1):NPART
        sampling_bar_p!(p,NBATHPART)
    end
end

##################################################

if INITFBATH == "GTE"
    const init_bathsampling! = sampling_GTE_bath!   # In states/GTE.jl file
elseif INITFBATH == "Polytrope"
    const init_bathsampling! = sampling_pol_bath!   # In states/Polytrope.jl file
elseif INITFBATH == "Plummer"
    const init_bathsampling! = sampling_plum_bath!  # In states/Plummer.jl file
elseif INITFBATH == "rand"
    const init_bathsampling! = sampling_rand_bath!
elseif INITFBATH == "bar"
    const init_bathsampling! = sampling_bar_bath!
end


##################################################
# Initialization of the test particles positions/velocities
# Caution ! There is an initialization before the warm up
# during which we assume only the bath interact with itself and
# all the "test" particles are massless during the warm up 
# EXCEPT for the "equal" mass sampling ! 
#
# Therefore in most cases, (x,v) initialization for test particles 
# is performed after warm up.
#
# In the "equal" mass sampling case, we consider that the test particles are some 
# bath particles (as they have the same mass) that the user want to follow.
# Therefore, using different DF samplings between the called bath and test particles 
# in this equal mass case is not recommended and (x,v) initialization for test particles 
# is performed before initialization.
##################################################

#####
# Sampling uniformly in a square of the phase space
# such that total energy is 3/4 * G * (MTOTBATH)^(2) * LAMBDA
# as for GTE (see Violent_relaxation MMA Notebook)
#####
function sampling_rand_test!()
    for p=1:NTESTPART
        sampling_rand_p!(p)
    end
end

#####
# Sampling uniformly between -4.0*Lambda and +4.0*Lambda in x
# with v = 0
#####
function sampling_unif_test!()
    for p=1:NTESTPART
        tabx[p] = 4.0*LAMBDA * 2.0*(rand() - 0.5)
        tabv[p] = v_from_meanDF_to_sim_RF(0.0)
    end
end


#####
# Sampling bar for phase mixing similar to wendy
# as for bath particles
#####
function sampling_bar_test!()
    for p=1:NTESTPART
        sampling_bar_p!(p,NTESTPART)
    end
end

#####
# Test particles are placed on one orbit in the phase space
# at uniformly distributed angle.
# Orbit : ra = lambda / 2
#####


function init_circle_test!()
    x0 = 4.0 / 3.0 * LAMBDA # nearly maximal diffusion at GTE in the Landau case
    E = pot(x0)

    tabx[1] = x_from_meanDF_to_sim_RF(x0)
    tabv[1] = v_from_meanDF_to_sim_RF(0.0)

    m = div(NTESTPART,2) # NTESTPART = 2*m + r
    r = NTESTPART%2 # NTESTPART = 2*m + r

    if m > 0
        for i = 1:(m-1)
            x = x0 * cos(2*pi*i/NTESTPART)
            potx = pot(x)
            v = sqrt(2*(E-potx))

            tabx[2*i] = x_from_meanDF_to_sim_RF(x)
            tabx[2*i+1] = x_from_meanDF_to_sim_RF(x)
            tabv[2*i] = v_from_meanDF_to_sim_RF(v)
            tabv[2*i+1] = v_from_meanDF_to_sim_RF(-v)
        end

        x = x0 * cos(2*pi*m/NTESTPART)
        potx = pot(x)
        v = sqrt(2*(E-potx))
        if r == 0
            tabx[2*m] = x_from_meanDF_to_sim_RF(x)
            tabv[2*m] = v_from_meanDF_to_sim_RF(v)
        else
            tabx[2*m] = x_from_meanDF_to_sim_RF(x)
            tabx[2*m+1] = x_from_meanDF_to_sim_RF(x)
            tabv[2*m] = v_from_meanDF_to_sim_RF(v)
            tabv[2*m+1] = v_from_meanDF_to_sim_RF(-v)
        end
    end
end

##################################################
if INITFTEST == "GTE"
    const testsampling! = sampling_GTE_test!    # In states/GTE.jl file
elseif INITFTEST == "Polytrope"
    const testsampling! = sampling_pol_test!    # In states/Polytrope.jl file
elseif INITFTEST == "Plummer"
    const testsampling! = sampling_plum_test!  # In states/Plummer.jl file
elseif INITFTEST == "rand"
    const testsampling! = sampling_rand_test!
elseif INITFTEST == "unif"
    const testsampling! = sampling_unif_test!
elseif INITFTEST == "bar"
    const testsampling! = sampling_bar_test!
elseif INITFTEST == "circle"
    const testsampling! = init_circle_test!
elseif INITFTEST == "oneplusGTE"
    const testsampling! = init_one_plus_GTE!    # In states/GTE.jl file
end

function donothing!()
end

if INITM == "equal"
    const init_testsampling_1! = testsampling!  # Before warm up
    const init_testsampling_2! = donothing!     # After warm up
else
    const init_testsampling_1! = donothing!     # Before warm up
    const init_testsampling_2! = testsampling!  # After warm up
end

##############################################################################################
# Initialization of the mass
# Caution ! There is an mass initialization before the warm up
# during which we assume only the bath interact with itself.
# All the "test" particles are massless during the warm up (and their position are 
# initialized after the warm up) EXCEPT for the "equal" mass sampling ! 
#
# In the "equal" mass sampling case, we consider that the test particles are some 
# bath particles (as they have the same mass) that the user want to follow.
# Therefore, using different DF samplings between the called bath and test particles 
# in this equal mass case is not recommended.
# 
# Therefore there is always 2 mass initialization functions (one called before the warm up 
# and one after)
#
# Extra assumption : the bath total mass is MTOTBATH (except for "equal" mass sampling). 
# Any massive test particle (e.g. for friction 
# illustration), add its mass to MTOTBATH for the total mass of the system MTOT[1]

# After any mass initialization, one has to update the value of the total mass in MTOT[1]
##############################################################################################

##################################################
# Equal mass particles
##################################################
function equal_mass_1!()
    for i=1:NPART
        tabm[i] = MTOTBATH/NPART # All particles are "bath" particles, 
    end
    MTOT[1]=sum(tabm)
end
function equal_mass_2!()
end

##################################################
# Equal mass particles
# + NTESTPART test particules without mass
##################################################
function equal_mass_plus_zero_1!() 
    for i=1:NTESTPART
        tabm[i] = 0.0
    end
    for i = (NTESTPART+1):NPART
        tabm[i] = MTOTBATH/NBATHPART
    end
    MTOT[1]=sum(tabm)
end
function equal_mass_plus_zero_2!() # No need to refix the test particle mass to 0
end

##################################################
# Equal mass particles
# + NTESTPART test particules much bigger
##################################################
function equal_mass_plus_big_1!() 
    for i=1:NTESTPART 
        tabm[i] = 0.0 # For warm up
    end
    for i = (NTESTPART+1):NPART # Mbath = MTOTBATH
        tabm[i] = MTOTBATH/NBATHPART
    end
    MTOT[1]=sum(tabm)
end
function equal_mass_plus_big_2!() 
    mtest =  10*MTOTBATH/NBATHPART 
    for i=1:NTESTPART
        tabm[i] = mtest # After warm up
    end
    MTOT[1]=sum(tabm)
end

##################################################
# 1 big + massless test particles (if more than 1 
# test particles)
# and equal mass bath
##################################################
function big_plus_nomass_1!() 
    for i=1:NTESTPART 
        tabm[i] = 0.0 # For warm up
    end
    for i = (NTESTPART+1):NPART # Mbath = MTOTBATH
        tabm[i] = MTOTBATH/NBATHPART
    end
    MTOT[1] = sum(tabm)
end
function big_plus_nomass_2!() 
    mtest =  sqrt(NBATHPART)*MTOTBATH/NBATHPART # mpart << sqrt(Nbath)*mpart << Nbath*mpart (Mbath) if N>>1
    tabm[1] = mtest
    # All other test particles are massless
    MTOT[1] = sum(tabm)
end


##################################################

if INITM == "equal"
    const init_mass_1! = equal_mass_1!
    const init_mass_2! = equal_mass_2!
elseif INITM == "testnomassp"
    const init_mass_1! = equal_mass_plus_zero_1!
    const init_mass_2! = equal_mass_plus_zero_2!
elseif INITM == "bigp"
    const init_mass_1! = equal_mass_plus_big_1!
    const init_mass_2! = equal_mass_plus_big_2!
elseif INITM == "bigandnomass"
    const init_mass_1! = big_plus_nomass_1!
    const init_mass_2! = big_plus_nomass_2!
end