##############################################################################################
# Usefull functions (potential, derivative ....)
# at Global Thermal Equilibrium (GTE)
##############################################################################################

##################################################
# On the same orbit (different position if there is more than one test particle)
# E = 0.5 * v^2 + psi(x) = cste
# Hypothesis : we are at GTE !!
# the potential psi(x) is then
# psi(x) = G * integral |x-x'| * rho(x') dx'
# with rho(x) = MTOTBATH / ( 2 * LAMBDA ) * sech^2( x / LAMBDA ).
# Then 
# psi(x) = MTOTBATH * G * LAMBDA * ln(2.0*cosh(x/LAMBDA)).
# And its derivative 
# psi'(x) = MTOTBATH * G * tanh(x/LAMBDA).
# We can invert this expression (taking only x>0) and find that
# x = LAMBDA * arccosh( 0.5 * exp( psi / (G*MTOTBATH*LAMBDA) ) ).
##################################################
function potential_GTE(x::Float64)
    return MTOTBATH * G * LAMBDA * log(2.0*cosh(x/LAMBDA))
end

function dpotdx_GTE(x::Float64)
    return MTOTBATH * G * tanh(x/LAMBDA)
end

function inv_pot_GTE(psi::Float64)
    return LAMBDA * acosh(0.5*exp(psi / (G*MTOTBATH*LAMBDA)))
end


##############################################################################################
# Initialization of the position and velocities
##############################################################################################


#####
# Sampling the Global Thermal Equilibrium (GTE)
#
# GTE distribution function
# f(x,v) = M/(2*sqrt(pi)*sigma*lambda)*exp(-v^2/sigma^2)*sech^2(x/lambda)
#        = M * g(x) * h(v)
# with 
#   g(x) = 1 / (2*lambda) * sech^2(x/lambda) 
# and  
#   h(v) = 1 / (sqrt(pi)*sigma)*exp(-v^2/sigma^2).
#
# To sample g (i.e. positions), we use its cumulative 
#   G(x) = int_{-\infty}^x g(x')dx' 
#        = 0.5 * ( 1 + tanh(x/lambda) ).
# Then, inversing the cumulative 
#   G^{-1}(u) = lambda * arctanh( 2 * u - 1 ),
# and sampling u uniformly, we got a simple sampling according to 
# the DF g.
#
# To sample h (i.e. velocities), we use randn() which follows 
# a normal centered distribution of std 1.
# As the normal (Gaussian) distribution N(m,s) is defined by
#   h_G (v) = 1 / ( sqrt(2*pi) * s ) * exp (-(v-m)^2/(2*s^2)),
# to simulate such a distribution, one just has to do
#       m + s * randn().
#
# In our case sigma is equivalent to sqrt(2) * s, (sigma = s * sqrt(2))
# then we sample
#       v = sigma / sqrt(2) * randn().
# 
# Finally, we convert the chosen velocity in the mean DF reference frame 
# into the simulation reference frame (which could be different, e.g. 
# when sampling test particle after a bath recentering).
#####
#####
# Sampling the particle p
#####
function sampling_GTE_p!(p::Int64)
    tabx[p] = LAMBDA * atanh( 2.0 * rand() - 1.0) # Sampling the sech^2(x/lambda) (cumulative : tanh)
    tabv[p] = v_from_meanDF_to_sim_RF( SIGMA * randn() / sqrt(2.0) )  # Sampling the e^(-v^2/sigma^2) (in the mean DF reference frame) and convert to simulation reference frame
end
#####
# Sampling the bath particles
#####
function sampling_GTE_bath!()
    for p=(NTESTPART+1):NPART
        sampling_GTE_p!(p)
    end
end

#####
# Sampling the Global Thermal Equilibrium (GTE)
# as for bath particles
#####
function sampling_GTE_test!()
    for p=1:NTESTPART
        sampling_GTE_p!(p)
    end
end


#####
# 1 test particle always placed at the same initial position
# + if more than 1 test particle : particles sampled as GTE
# Useful for one massive first test particles and 
# massless others
#####

function init_one_plus_GTE!()
    x0 = 4.0 / 3.0 * LAMBDA # nearly maximal diffusion at GTE in the Landau case
    E = pot(x0)

    tabx[1] = x0
    tabv[1] = v_from_meanDF_to_sim_RF(0.0) # Coordinate w.r.t. the simulation reference frame

    for p = 2:NTESTPART # May be void
        sampling_GTE_p!(p)
    end
end