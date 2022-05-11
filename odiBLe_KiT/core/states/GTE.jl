############################################################################
#
# SPECIFIC TO GTE
#
############################################################################


##################################################
# Mean potential for the GTE
##################################################
function psi_GTE(ra::Float64)
    return Mtot * G * Lambda * log(2.0 * cosh(ra/Lambda))
end
##########
# Derivative of the mean potential for the GTE
##########
function dpsidra_GTE(ra::Float64)
    return Mtot * G * tanh(ra/Lambda)
end
##########
# Second derivative of the mean potential for the GTE
##########
function d2psidra2_GTE(ra::Float64)
    return Mtot * G * (sech(ra/Lambda))^(2) / Lambda 
end
##########
# Third derivative of the mean potential for the GTE
##########
function d3psidra3_GTE(ra::Float64)
    return - 2.0 * Mtot * G / (Lambda)^(2) * tanh(ra/Lambda) * (sech(ra/Lambda))^(2)
end
##########
# Fourth derivative of the mean potential for the GTE
##########
function d4psidra4_GTE(ra::Float64)
    return - 2.0 * Mtot * G / (Lambda)^(3) * (sech(ra/Lambda))^(2) * (1 - 3.0 * (tanh(ra/Lambda))^(2))
end
##########
# DF for the GTE, as a function of E
##########
function DFfromE_GTE(E::Float64)
    return 2.0 * Mtot * exp(- 2.0 * E / (Mtot * G * Lambda)) / (sqrt(pi) * Sigma * Lambda)
end
##########
# dDF/dE for the GTE, as a function of E
##########
function dDFdEfromE_GTE(E::Float64)
    return - 4.0 * exp(- 2.0 * E / (Mtot * G * Lambda)) / (sqrt(pi) * G * Sigma * (Lambda)^(2) )
end


# In all these functions, ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
const psi = psi_GTE                 # Mean potential
const dpsidra = dpsidra_GTE         # Mean potential derivative
const d2psidra2 = d2psidra2_GTE     # Mean potential second derivative
const d3psidra3 = d3psidra3_GTE     # Mean potential third derivative
const d4psidra4 = d4psidra4_GTE     # Mean potential fourth derivative
const DFfromE = DFfromE_GTE         # DF as a function of E
const dDFdEfromE = dDFdEfromE_GTE   # dDF/dE as a function of E