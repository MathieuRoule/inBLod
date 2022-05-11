############################################################################
#
# Potential and other functions specific to the considered QSS
#
############################################################################
include("./states/"*STATE*".jl")


############################################################################
#
# MAPPING FUNCTION x(u) = ra * f(u)
#
############################################################################
# The mapping function must satisfy :
# +  f : [-1,1] --> [-1,1] odd function
# +  f(+1) = +1 (and f(-1) = -1)
# +  f'(1) = f'(-1) = 0 
############################################################################
include("./mapping/"*MAPPING*".jl")


############################################################################
#
# GENERIC FUNCTIONS
#
############################################################################
#####
# Value for change between true expression and taylor expansion
# both for the anomaly u -> 1 and the apocentre ra -> 0
#####
const eps_u = 10^(-3)
const eps_ra = 10^(-3)
##################################################
# Very usefull function
# df_sqrtpsi(ra,u) = f'(u) / sqrt( psi(ra) - psi(ra*f(u)) )
# ATTENTION, DL in u = +- 1
##################################################
function df_sqrtpsi(ra::Float64,u::Float64) 
    @assert (abs(u) <= 1.0+10^(-14)) "abs(u) greater than 1 in df_sqrtpsi"
    if (abs(1-abs(u)) < eps_u)
        if u > 0.0
            return df_sqrtpsi_bounds(ra,u,true)
        else
            return df_sqrtpsi_bounds(ra,u,false)
        end
    else
        return dfdu(u) / sqrt(psi(ra) - psi(ra * fu(u)))
    end
end

##################################################
# DL (Taylor's polynomial approx.) of
# df_sqrtpsi(ra,u) = f'(u) / sqrt( psi(ra) - psi(ra*f(u)) )
# in u = +- 1
# at second order
# if l == False : DL in -1 
# if l == True : DL in +1
##################################################
function df_sqrtpsi_bounds(ra::Float64,u::Float64,l::Bool)
    ul = l ? 1.0 : -1.0
    pref = - ul
    xl = ra * ul

    f2 = d2fdu2(ul)
    f3 = d3fdu3(ul)
    f4 = d4fdu4(ul)

    p1 = dpsidra(xl)
    p2 = d2psidra2(xl)

    den = sqrt(-ra * p1 * f2)

    a = sqrt(2.0) * f2
    b = sqrt(2.0) * f3 / 3.0
    c = (3.0 * (p1 * f2 * f4 - ra * p2 * (f2)^(3)) - p1 * (f3)^(2)) / (12.0 * sqrt(2.0) * p1 * f2)

    return pref / den * ( a + b * (u - ul) + c * (u - ul)^(2) )
end

##################################################
# Integration using RK4 scheme
# When integrate for u between 0 and 1.
##################################################
function integrate_u_01(fun,ra::Float64)
    Delta_u = 1.0 / K_MF
    u = 0.0
    res = Delta_u * fun(ra,u) / 6.0

    for k = 1:(K_MF-1)
        u += 0.5*Delta_u
        res += 2.0 * Delta_u * fun(ra,u) / 3.0
        u += 0.5*Delta_u
        res += 1.0 * Delta_u * fun(ra,u) / 3.0
    end
    u += 0.5*Delta_u
    res += 2.0 * Delta_u * fun(ra,u) / 3.0
    u += 0.5*Delta_u
    res += 1.0 * Delta_u * fun(ra,u) / 6.0

    return res
end

##################################################
# Action from the apocentre
# WARNING : Could be false at high ra
##################################################
function dftimessqrtpsi(ra::Float64,u::Float64)
    @assert (abs(u) <= 1.0+10^(-14)) "abs(u) greater than 1 in dftimessqrtpsi"
    return dfdu(u) * sqrt(psi(ra) - psi(ra * fu(u)))
end
function Jfromra(ra::Float64)
    if eps_ra < ra 
        return 2.0 * sqrt(2.0) * ra * integrate_u_01(dftimessqrtpsi,ra) / pi
    else
        return 0.5 * sqrt(d2psidra2(0.0)) * (ra)^(2) 
    end
end

##################################################
# Frequency
# as a function of the orbit's apocentre
# WARNING : Could be false at high ra
##################################################
function Omegafromra(ra::Float64)
    if eps_ra < ra 
        inv_omega = sqrt(2.0) * ra * integrate_u_01(df_sqrtpsi,ra) / pi
        return 1.0 / inv_omega
    else
        return sqrt(d2psidra2(0.0)) * (1 +  d4psidra4(0.0) * (ra)^(2) / (16.0 * d2psidra2(0.0)) ) 
    end
end


##################################################
# Function to integrate for u : 0 -> 1 to get dOmega/dra
# (up to a prefactor)
##################################################
function fundodr(ra::Float64,u::Float64) # Hypothesis u >= 0 (DL only in u = +1)
    @assert (abs(u) <= 1.0+10^(-14)) "abs(u) greater than 1 in fundodr"
    if (abs(1-abs(u)) < eps_u) 
        return fundodr_bounds(ra,u) 
    else
        return dfdu(u) * (dpsidra(ra) - fu(u) * dpsidra(ra * fu(u))) / (sqrt(psi(ra) - psi(ra * fu(u))))^(3)
    end
end

##################################################
# DL (Taylor's polynomial approx.) of 
# the function to integrate to get dOmega/dra
# in u = +1
##################################################
function fundodr_bounds(ra::Float64,u::Float64)
    f2 = d2fdu2(1.0)
    f3 = d3fdu3(1.0)

    p1 = dpsidra(ra)
    p2 = d2psidra2(ra)

    pref = sqrt(2.0) * (p1 + ra * p2) / (sqrt(-ra * p1 * f2))^(3)
    
    a = (f2)^(2)
    b = f2 * f3 / 3.0

    return pref * ( a + b * (u - 1.0) )
end

##################################################
# Value of dOmega/dra
##################################################
function dOmegadra(ra::Float64)
    Omega = Omegafromra(ra)
    if eps_ra < ra
        return - Omega / ra + ra * (Omega)^(2) * integrate_u_01(fundodr,ra) / (sqrt(2.0) * pi)
    else
        return d4psidra4(0.0) * ra / (8.0 * sqrt(d2psidra2(0.0)))          
    end
end


##################################################
# Energy from the apocentre
##################################################
function Efromra(ra::Float64)
    return psi(ra) # Orbital energy
end
##################################################
# Jacobian |dJ/dra|
# |dJ/dra| = |dJ/dE| * |dE/dra| = |psi'(ra) / Omega(ra)|
##################################################
function JacJtora(ra::Float64)
    return abs(dpsidra(ra) / Omegafromra(ra))
end

##################################################
# Jacobian dtheta/du for a given orbit
##################################################
function dthetadu(ra::Float64,u::Float64)
    return Omegafromra(ra) * ra * df_sqrtpsi(ra,u) / sqrt(2.0) 
end
##################################################
# Jacobian dtheta/du for a given orbit
# Divided by Omega (so we don't call Omegafromra too 
# many times)
##################################################
function dthetadu_omega(ra::Float64,u::Float64)
    return ra * df_sqrtpsi(ra,u) / sqrt(2.0) 
end

##################################################
# DF from ra
# ATTENTION, normalised so that INT[DF,{x,v}] = Mtot
##################################################
function DF(ra::Float64)
    E = Efromra(ra) # Orbital energy
    return DFfromE(E)
end
##################################################
# dDF/dJ from ra
# We use the relation dE/dJ=Omega
##################################################
function dDFdJ(ra::Float64)
    E = Efromra(ra) # Orbital energy
    Omega = Omegafromra(ra) # Orbital frequency
    return dDFdEfromE(E)*Omega
end
