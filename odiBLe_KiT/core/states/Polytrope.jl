############################################################################
#
# SPECIFIC TO POLYTROPE QSS
#
############################################################################


##################################################
# Mean potential for the Polytrope
##################################################
function psi_poly(ra::Float64)
    if (abs(ra)<Alpha)
        return G * Mtot * Alpha * ( 1.0 - 2.0/pi * cos((pi*ra)/(2.0*Alpha)) )
    else
        return G * Mtot * abs(ra)
    end 
end
##########
# Derivative of the mean potential for the Polytrope
##########
function dpsidra_poly(ra::Float64)
    if (abs(ra)<Alpha)
        return G * Mtot * sin((pi*ra)/(2.0*Alpha))
    else
        return G * Mtot * sign(ra)
    end 
end
##########
# Second derivative of the mean potential for the Polytrope
##########
function d2psidra2_poly(ra::Float64)
    if (abs(ra)<Alpha)
        return G * Mtot * pi / (2.0*Alpha) * cos((pi*ra)/(2.0*Alpha))
    else
        return 0.0
    end 
end
##########
# Third derivative of the mean potential for the Polytrope
##########
function d3psidra3_poly(ra::Float64)
    if (abs(ra)<Alpha)
        return - G * Mtot * (pi / (2.0*Alpha))^(2.0) * sin((pi*ra)/(2.0*Alpha))
    else
        return 0.0
    end 
end
##########
# Fourth derivative of the mean potential for the Polytrope
##########
function d4psidra4_poly(ra::Float64)
    if (abs(ra)<Alpha)
        return - G * Mtot * (pi / (2.0*Alpha))^(3.0) * cos((pi*ra)/(2.0*Alpha))
    else
        return 0.0
    end 
end
##########
# Energy which are populated in Polytrope
##########
function ispopulated_poly(E::Float64)
    if 1.0 - 2.0/pi <= E / (G * Mtot * Alpha) < 1.0
        return true
    else
        return false
    end
end
##########
# DF for the Polytrope, as a function of E
##########
function DFfromE_poly(E::Float64)
    if ispopulated_poly(E)
        Mtot * pi / (8.0 * Alpha) * sqrt(2.0 / (G * Mtot * Alpha) * (1.0 - E /(G * Mtot * Alpha)))
    else
        return 0.0
    end
end
##########
# dDF/dE for the Polytrope, as a function of E
##########
function dDFdEfromE_poly(E::Float64)
    if ispopulated_poly(E)
        - pi / (8.0 * G * (Alpha)^(2.0) * sqrt(2.0 * G * Mtot * Alpha)) * 1.0 / sqrt(1.0 - E /(G * Mtot * Alpha))
    else
        return 0.0
    end
end


# In all these functions, ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
const psi = psi_poly                # Mean potential
const dpsidra = dpsidra_poly        # Mean potential derivative
const d2psidra2 = d2psidra2_poly    # Mean potential second derivative
const d3psidra3 = d3psidra3_poly    # Mean potential third derivative
const d4psidra4 = d4psidra4_poly    # Mean potential fourth derivative
const DFfromE = DFfromE_poly        # DF as a function of E
const dDFdEfromE = dDFdEfromE_poly  # dDF/dE as a function of E