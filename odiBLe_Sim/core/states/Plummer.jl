##############################################################################################
# Usefull functions (potential, derivative ....)
# for the Plummer QSS (density is a Lorentzian to the power 3/2)
##############################################################################################

#####
# Length scale parameter such that Plummer state has 
# the same energy as the GTE 
####
ALPHAPLUMMER = 2.0 / pi * LAMBDA 

##################################################
# On the same orbit (different position if there is more than one test particle)
# E = 0.5 * v^2 + psi(x) = cste
#
# QUASI-STATIONNARY STATE --> Plummer : Lorentzian to the power 3/2
#
# the potential psi(x) is then
# psi(x) = G * integral |x-x'| * rho(x') dx'
# with rho(x) = MTOTBATH / ( 2 * ALPHAPLUMMER ) * 1  /  ( 1 + (x/ALPHAPLUMMER)^2 )^(3/2),
# 
# Then 
# psi(x) = G * MTOTBATH * ALPHAPLUMMER * sqrt( 1 + (x/ALPHAPLUMMER)^2 )
# with ALPHAPLUMMER = 2 / pi * LAMBDA so this QSS has the same total energy and mass as 
# the corresponding LAMBDA GTE state.
#
# And its derivative
# psi'(x) = G * MTOTBATH * x/ALPHAPLUMMER * 1 / sqrt( 1 + (x/ALPHAPLUMMER)^2 ).
#
# We can invert the potential expression and find that
# x = ALPHAPLUMMER * sqrt( (psi / G * MTOTBATH * ALPHAPLUMMER)^2 - 1 ),
# for psi >= G*MTOTBATH*ALPHAPLUMMER
##################################################
function potential_plum(x::Float64)
    return G * MTOTBATH * ALPHAPLUMMER * sqrt( 1.0 + (x/ALPHAPLUMMER)^2 )
end

function dpotdx_plum(x::Float64)
    return G * MTOTBATH * x/ALPHAPLUMMER * 1.0 / sqrt( 1.0 + (x/ALPHAPLUMMER)^2 )
end

function inv_pot_plum(psi::Float64)
    if psi < G * MTOTBATH * ALPHAPLUMMER
        return NaN
    else 
        return ALPHAPLUMMER * sqrt( (psi / G * MTOTBATH * ALPHAPLUMMER)^2 - 1 )
    end
end


##############################################################################################
# Initialization of the position and velocities
##############################################################################################

##################################################
# Sampling the Plummer Quasi-Stationnary State
#
# Plummer distribution function
# f(x,v) = !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# To sample this function we first sample the positions according to the density 
# rho(x) = 1 / ( 2 * ALPHAPLUMMER ) * 1  /  ( 1 + (x/ALPHAPLUMMER)^2 )^(3/2) (*MTOTBATH -> normalized at 1 for the sampling)
#
# The cumulative of this density is 
# cumrho(x) = 1/2  *  (  1  +  (x/ALPHAPLUMMER) / sqrt(1 + (x/ALPHAPLUMMER)^2)  )  (*MTOTBATH idem)
# Then the inverse of its cumulative is 
# invcumrho(u) = ALPHAPLUMMER *  (2 * u - 1) / sqrt( 1 - (2*u-1)^2 )
#
# We then have a simple sampling for the positions, using a uniform distribution and the inverse cumulative.
#
# For the velocities, they have to be linked to the position of the particle (as the Hamiltonian is not separable as 
# for the GTE).
# 
# We use rejection sampling with the probability distribution of y = sqrt( gamma/beta(x) ) * v 
# with 
#       gamma = 1/(2*G*MTOTBATH*ALPHAPLUMMER)
# and
#       beta(x) = sqrt( 1 + (x/ALPHAPLUMMER)^2 ),
#
# f(y) = 15/16 * ( 1 + y^2 )^(-7/2).
#
# To sample this PDF we dominate it with 
# 
# g(y) =  |y|^(-7)    if |y| > yc,  
#      =   yc^(-7)    if not
# with yc = (7/3)^(1/6) such that g is normalized.
#
# With this PDF g, one has
#       f <= C0 * g      with C0 = (35*yc)/16 < 2.52
#
# We sample g using its cumulative 
#
# cumg(y) =     y^(-6)/6    if  y < -yc
#         = 1 - y^(-6)/6    if  y >  yc
#         = 1/2 + yc^(-7) * y else
#
# which inverse is 
#
# invcumg(u)  =  -(6 u)^(-1/6)      if u < 1/14
#             =  (6*(1-u))^(-1/6)   if u > 13/14
#             =  (u-1/2) * yc^7     else
#
# Finally, we convert the chosen velocity in the mean DF reference frame 
# into the simulation reference frame (which could be different, e.g. 
# when sampling test particle after a bath recentering).
##################################################

const gammaplum = 1.0/(2.0*G*MTOTBATH*ALPHAPLUMMER) # Characteristic velocity squared
#####
# Function that control the spread of velocity distribution knowing the position x
#####
function betaplum(x::Float64)
    return sqrt( 1 + (x/ALPHAPLUMMER)^2 )
end
#####
# PDF to sample depending on y = sqrt( gamma/beta(x) ) * v 
#####
function pdfyplum(y::Float64)
    return 15.0/16.0 * ( sqrt(1 + (y)^(2)) )^(-7)
end
#####
# Dominating PDF for rejection method (to multiply by a constant)
#####
const ycgdomplum = (7.0/3.0)^(1.0/6.0) # Cut between power law and flat in the center
function pdfgdomplum(y::Float64)
    if abs(y) < ycgdomplum
        return (ycgdomplum)^(-7)
    else
        return (abs(y))^(-7)
    end
end
const Cdomplum = 2.52 # Constant such as f < C * g
#####
# Cumulative of the dominating function
#####
function cumGdomplum(y::Float64)
    if y < - ycgdomplum
        return (y)^(-6) / 6.0
    elseif y < ycgdomplum
        return 0.5 + (ycgdomplum)^(-7) * y
    else
        return 1.0 - (y)^(-6) / 6.0
    end
end
#####
# Inverse of the cumulative of the dominating function
#####
function invcumGdomplum(u::Float64)
    if (u <= 0.0) || (u >= 1.0)
        return 0.0
    elseif u < 1.0/14.0
        return -(6.0*u)^(-1.0/6.0)
    elseif u > 13.0/14.0
        return (6.0*(1.0 - u))^(-1.0/6.0)
    else
        return (u - 0.5)*(ycgdomplum)^(7)
    end
end

#####
# Sampling the particle p
#####
function sampling_plum_p!(p::Int64)
    # Sampling the position
    u = rand()
    x = ALPHAPLUMMER *  (2.0*u - 1.0) / sqrt( 1.0 - (2.0*u - 1.0)^(2) ) # Sampling the density (cumulative : y/sqrt(1-y^2))
    tabx[p] = x

    # Sampling the velocity by rejection 
    reject = true
    while reject
        yprop = invcumGdomplum(rand())
        utest = rand()
        if utest < pdfyplum(yprop) / (Cdomplum * pdfgdomplum(yprop))
            # Converting the y proposal into its associated velocity v = sqrt( beta(x) / gamma ) * y
            tabv[p] = v_from_meanDF_to_sim_RF( sqrt( betaplum(tabx[p]) / gammaplum) * yprop ) # Coordinate w.r.t. the simulation reference frame
            reject = false
        end
    end
end
#####
# Sampling the bath particles
#####
function sampling_plum_bath!()
    for p=(NTESTPART+1):NPART
        sampling_plum_p!(p)
    end
end

#####
# Sampling the test particles
#####
function sampling_plum_test!()
    for p=1:NTESTPART
        sampling_plum_p!(p)
    end
end
