############################################################################
#
# SPECIFIC TO 1D PLUMMER QSS
#
############################################################################
# which correspond to a lorentzian function to the power 3/2
# rho(x)  =  M / (2*alpha) * (1 + (x/alpha)^2)^(-3/2)
############################################################################

##################################################
# Mean potential for the Plummer state
##################################################
function psi_Plummer(ra::Float64)
    return G * Mtot * Alpha * sqrt( 1.0 + (ra/Alpha)^(2) )
end
##########
# Derivative of the mean potential for the Plummer state
##########
function dpsidra_Plummer(ra::Float64)
    return G * Mtot * (ra/Alpha) * 1.0 / sqrt( 1.0 + (ra/Alpha)^(2) )
end
##########
# Second derivative of the mean potential for the Plummer state
##########
function d2psidra2_Plummer(ra::Float64)
    return G * Mtot / Alpha * 1.0 / ( (1.0 + (ra/Alpha)^(2)) * sqrt(1.0 + (ra/Alpha)^(2)) )
end
##########
# Third derivative of the mean potential for the Plummer state
##########
function d3psidra3_Plummer(ra::Float64)
    return - 3.0 * G * Mtot / (Alpha)^(2) * ra/Alpha * 1.0 / ( (1.0 + (ra/Alpha)^(2))^(2) * sqrt(1.0 + (ra/Alpha)^(2)) )
end
##########
# Fourth derivative of the mean potential for the Plummer state
##########
function d4psidra4_Plummer(ra::Float64)
    return - 3.0 * G * Mtot / (Alpha)^(3) * (1.0 - 4.0*(ra/Alpha)^(2)) * 1.0 / ( (1.0 + (ra/Alpha)^(2))^(3) * sqrt(1.0 + (ra/Alpha)^(2)) )
end
##########
# DF for the Plummer state, as a function of E
##########
function DFfromE_Plummer(E::Float64)
    return ( 15.0 / (32.0 * sqrt(2.0)) ) * (G)^(3) * (Mtot)^(4) * (Alpha)^(2) / ( (E)^(3) * sqrt(E) ) 
end
##########
# dDF/dE for the Plummer state, as a function of E
##########
function dDFdEfromE_Plummer(E::Float64)
    return - ( 105.0 / (64.0 * sqrt(2.0)) ) * (G)^(3) * (Mtot)^(4) * (Alpha)^(2) / ( (E)^(4) * sqrt(E) ) 
end


# In all these functions, ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
const psi = psi_Plummer                 # Mean potential
const dpsidra = dpsidra_Plummer         # Mean potential derivative
const d2psidra2 = d2psidra2_Plummer     # Mean potential second derivative
const d3psidra3 = d3psidra3_Plummer     # Mean potential third derivative
const d4psidra4 = d4psidra4_Plummer     # Mean potential fourth derivative
const DFfromE = DFfromE_Plummer         # DF as a function of E
const dDFdEfromE = dDFdEfromE_Plummer   # dDF/dE as a function of E
