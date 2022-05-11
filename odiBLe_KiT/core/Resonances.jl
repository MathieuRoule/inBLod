##################################################
# Defining the functions used for finding zeros
##################################################
# Function used to find zero
# It is a simple bisection algorithm, but it makes no allocations and is sufficiently fast
# It allows us not to have to use the Roots library that makes a lot of allocations
# The optional tolerances are set to the same as the ones I found in Roots.jl.
# @ATTENTION, the tolerances are most likely overkill -- it may prevent convergence for high-order resonances
# @IMPROVE -- it could be a good idea to put a counter of iterations, to prevent the algorithm from getting stuck?
##################################################
function bisection(fun, xl::Float64, xu::Float64, tolx::Float64=1000.0*eps(Float64), tolf::Float64=1000.0*eps(Float64))
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET"
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if (abs(xu-xl) <= tolx) # The considered bracket is smaller than the tolerance
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
end


##################################################
# Resonances research 
##################################################
#####
# Bounds used for the bracketing of the bisection
#####
# Bounds in ra
const RA_RESONANCE_MIN = 0.0
#####
# Bounds in Omega
# ATTENTION, ra->Omega(ra) is a decreasing function
const OMEGA_RESONANCE_MAX = Omegafromra(RA_RESONANCE_MIN)


################################################
# Maximal apocentre to look for a resonance for GTE State
#
# Using the fact that the orbital time of a particle
# will always be lower than if the total mass was
# concentrated in 0, we can deduce a function
# Om_+(ra) which dominate Omegafromra(ra) whatever 
# the distribution function is.
# 2 * pi / Om_+(ra) = 2 * sqrt(2) * int_0^ra dx / sqrt(G*M*ra - G*M*|x|)
#                   = 2 * sqrt(2) / sqrt(G*M*ra) * int_0^1 radu / sqrt(1-u)
#                   = (2 * sqrt(2) * sqrt(ra) / sqrt(G*M)) * 2.
# Then, we have
# Om_+(ra) = pi * sqrt(G*M) / (2 * sqrt(2) * sqrt(ra)) >= Omegafromra(ra)
# and Omegafromra(ra) and Om_+(ra) are decreasing functions of ra,
# then for a given omega, the true corresponding ra is lower than
# ra_sup such that Om_+(ra_sup) = omega.
# Therefore, a good superior value for ra is 
# ra_sup = G * M * pi^2 / (8 * omega^2)
################################################
function get_ramax_GTE(omega::Float64)
    return G * Mtot * (pi)^2 / (8.0 * (omega)^(2))
end
################################################
# Maximal apocentre to look for a resonance for Polytrope n=1 State
################################################
function get_ramax_pol(omega::Float64)
    return Alpha
end
################################################
# Maximal apocentre to look for a resonance for Plummer state
# Same as GTE with the same analysis
################################################
const get_ramax_Plummer = get_ramax_GTE
#####
if STATE == "GTE"
    const get_ramax = get_ramax_GTE
elseif STATE == "Polytrope"
    const get_ramax = get_ramax_pol 
elseif STATE == "Plummer"
    const get_ramax = get_ramax_Plummer  
end

##################################################
# For a given value of omega and np
# returns the value of rap, if it exists,
# so that kp*Omega(rap) = omega
##################################################
function get_Resonance(omega::Float64,kp::Int64)
    # The resonance condition can be solved
    ra_max = get_ramax(omega/(kp))
    omega_min = Omegafromra(ra_max)
    if (omega_min < omega/(kp) < OMEGA_RESONANCE_MAX)
        return bisection(rap -> kp*Omegafromra(rap)-omega, RA_RESONANCE_MIN, ra_max)
    # The resonance condition cannot be solved
    else
        return nothing
    end
end



########################################################################
# Count resonances
# Count the number of resonances which will be explored when computing
# the diffusion or friction coefficients at a given apocenter 
########################################################################
function count_resonances(ra::Float64) 
    count = 0 # Counter 
    #####
    val_Omega = Omegafromra(ra) # Current value of the resonance frequency
    #####
    for k=1:NMAX # Loop over the k. ATTENTION, we limit ourselves to k>0
        kdotomega = k*val_Omega # Current resonance frequency
        #####
        # Loop over the kp
        # ATTENTION, we limit ourselves to kp>0
        # ATTENTION, we limit ourselves to (k-kp) even
        # @IMPROVE -- the starting point for kp could be written in a more readable form
        for kp=(mod(k+1,2)+1):2:NMAX # Loop over np. ATTENTION, the `step' is 2
            rap = get_Resonance(kdotomega,kp) # Solving for the resonance condition
            #####
            if (rap != nothing) # We have found a solution to the resonance condition
                count += 1
            end
            #####
        end
    end
    #####
    return count # Output
end
