##############################################################################################
# Usefull functions (potential, derivative ....)
# for the Polytrope 1D n=1 QSS
##############################################################################################

#####
# Length scale parameter such that Lorentz32 state has 
# the same energy as the GTE 
####
ALPHAPOLY = 2.0 * LAMBDA 

##################################################
# On the same orbit (different position if there is more than one test particle)
# E = 0.5 * v^2 + psi(x) = cste
#
# QUASI-STATIONNARY STATE --> Polytrope 1D n=1
#
# the potential psi(x) is then
# psi(x) = G * integral |x-x'| * rho(x') dx'
# with rho(x) = MTOTBATH * pi / ( 4 * ALPHAPOLY ) * cos( pi / 2 * x / ALPHAPOLY ),
# if -ALPHAPOLY =< x =< ALPHAPOLY, 0 elsewhere.
# 
# Then 
# psi(x) = G * MTOTBATH * ALPHAPOLY * [1 - 2/pi * cos( pi / 2 * x / ALPHAPOLY )]
# if -ALPHAPOLY =< x =< ALPHAPOLY, G * MTOTBATH * |x| elsewhere,
# with ALPHAPOLY = 2 * LAMBDA so this QSS has the same total energy and mass as 
# the corresponding LAMBDA GTE state.
#
# And its derivative
# psi'(x) = G * MTOTBATH * sin( pi / 2 * x / ALPHAPOLY )
# if -ALPHAPOLY =< x =< ALPHAPOLY, G * MTOTBATH * sign(x) elsewhere
#
# We can invert the potential expression and find that
# x = ALPHAPOLY * 2/pi * arccos( pi/2 * [1 - psi/(G*MTOTBATH*ALPHAPOLY)]),
# if psi =< G*MTOTBATH*ALPHAPOLY, psi / (G*MTOTBATH) elsewhere.
##################################################
function potential_pol(x::Float64)
    if abs(x) > ALPHAPOLY
        return G * MTOTBATH * abs(x)
    else
        return G * MTOTBATH * ALPHAPOLY * ( 1.0 - 2.0/pi * cos( (pi*x)/(2.0*ALPHAPOLY) ) )
    end
end

function dpotdx_pol(x::Float64)
    if x > ALPHAPOLY
        return G * MTOTBATH
    elseif x < - ALPHAPOLY
        return - G * MTOTBATH
    else
        return G * MTOTBATH * sin( (pi*x)/(2.0*ALPHAPOLY) )
    end
end

function inv_pot_pol(psi::Float64)
    if psi < (1.0 - 2.0/pi) * G * MTOTBATH * ALPHAPOLY
        return NaN
    elseif psi > G * MTOTBATH * ALPHAPOLY
        return psi / (G * MTOTBATH)
    else 
        return ALPHAPOLY * 2.0 / pi * acos( pi/2.0 * (1.0 - psi/(G*MTOTBATH*ALPHAPOLY) ) )
    end
end

##############################################################################################
# Initialization of the position and velocities
##############################################################################################

##################################################
# Sampling the Polytrope Quasi-Stationnary State
#
# Polytrope distribution function
# f(x,v) = 1/sqrt[2*G*MTOTBATH*ALPHAPOLY] * sqrt[ 2/pi*cos(pi/2*x/ALPHAPOLY) - v^2/(2*G*M*ALPHAPOLY) ] / cos(pi/2*x/ALPHAPOLY)
# if abs(x) =< ALPHAPOLY and abs(v) =< sqrt[ 4*G*MTOTBATH*ALPHAPOLY/pi * cos(pi/2*x/ALPHAPOLY) ], 0 elsewhere.
#
# To sample this function we first sample the positions according to the density 
# rho(x) = pi/(4*ALPHAPOLY) * cos(pi/2*x/ALPHAPOLY) (*MTOTBATH -> normalized at 1 for the sampling)
# if abs(x) =< ALPHAPOLY, 0 elsewhere.
# The cumulative of this density is 
# cumrho(x) = 0.5 * [ 1 + sin( (pi*x)/(2*ALPHAPOLY) ) ]  (*MTOTBATH idem)
# Then the inverse of its cumulative is 
# invcumrho(u) = ALPHAPOLY * 2/pi *arcsin(2*u - 1)
#
# We then have a simple sampling for the positions, using a uniform distribution and the inverse cumulative.
#
# For the velocities, they have to be linked to the position of the particle (as the Hamiltonian is not separable as 
# for the GTE).
# 
# We use rejection sampling with the probability distribution of v given x 
# (we here suppose that x is between -ALPHAPOLY and ALPHAPOLY - which would be the case)
# 
# P(v|x) = 1/sqrt(2*G*M*ALPHAPOLY) * sqrt[ 2/pi * cos(pi/2*x/ALPHAPOLY) - v^2/(2*G*M*ALPHAPOLY) ] / cos(pi/2*x/ALPHAPOLY)
# if abs(v) =< sqrt[ 4*G*M*ALPHAPOLY/pi * cos(pi/2*x/ALPHAPOLY) ] = vmax(x)  , 0 elsewhere
# 
# We use a uniform sampling on [-vmax(x), vmax(x)[ and the constant (depending on x not on v)
# C = P(0|x) = 1/sqrt(pi*G*M*ALPHAPOLY) / sqrt[ cos(pi/2*x/ALPHAPOLY) ] 
# to dominate the true distribution function
# as P(v|x) =< P(0|x) for all x and v.
#
# Finally, we convert the chosen velocity in the mean DF reference frame 
# into the simulation reference frame (which could be different, e.g. 
# when sampling test particle after a bath recentering).
##################################################
#####
# Maximal velocity for a given position x in the polytrope sampling
#####
function pol_vmax(x::Float64)
    if abs(x) >= ALPHAPOLY
        return 0.0
    else 
        return sqrt( 4.0 * G*MTOTBATH*ALPHAPOLY / pi * cos( pi/2.0 * x/ALPHAPOLY ) )
    end
end
#####
# Conditional probability of v given x for the polytrope sampling
#####
function pol_Pvgivenx(x::Float64,v::Float64)
    if v >= pol_vmax(x)
        return 0.0
    else
        cx = cos(pi/2.0 * x/ALPHAPOLY)
        return 1.0 / sqrt(2.0*G*MTOTBATH*ALPHAPOLY) * sqrt( 2.0/pi * cx - v^2 / (2.0*G*MTOTBATH*ALPHAPOLY) ) / cx
    end
end
#####
# Sampling the particle p
#####
function sampling_pol_p!(p::Int64)
    # Position
    x = ALPHAPOLY * 2.0/pi * asin( 2.0*rand() - 1.0 ) # Sampling the density (cumulative : arcsin)
    tabx[p] = x

    # Velocity
    vmax = pol_vmax(x)
    if vmax == 0.0 
        error("Warning : Sampled velocity forced to 0")
        tabv[p] = v_from_meanDF_to_sim_RF(0.0)
        reject=false
    else
        reject=true
    end 
    Cmaj = pol_Pvgivenx(x,0.0)      # Majorant to the distribution function
    while reject
        v =  vmax * ( 2.0 * rand() - 1.0 )  # take u uniformly between -vmax and +vmax
        ut = rand()                         # take ut uniformly between 0 and 1
        if ut <= pol_Pvgivenx(x,v)/ Cmaj     # We accept the velocity u with probability P(v|x)
            tabv[p] = v_from_meanDF_to_sim_RF(v) # Coordinate w.r.t. the simulation reference frame
            reject=false
        end
    end
end
#####
# Sampling the bath particles
#####
function sampling_pol_bath!()
    for p=(NTESTPART+1):NPART
        sampling_pol_p!(p)
    end
end

#####
# Sampling the QSS Polytrope 1D n=1
# as for bath particles
#####
function sampling_pol_test!()
    for p=1:NTESTPART
        sampling_pol_p!(p)
    end
end