##################################################
# Polynomial Effective anomaly
##################################################
function fu_pol(u::Float64)
    return u*(1.5 - 0.5*u^(2)) # Output
end
##########
# Derivative of the polynomial effective anomaly
##########
function dfdu_pol(u::Float64)
    return 1.5*(1.0 - u^(2)) # Output
end
##########
# Second derivative of the polynomial effective anomaly
##########
function d2fdu2_pol(u::Float64)
    return - 3.0 * u # Output
end
##########
# Third derivative of the polynomial effective anomaly
##########
function d3fdu3_pol(u::Float64)
    return - 3.0 # Output
end
##########
# Fourth derivative of the polynomial effective anomaly
##########
function d4fdu4_pol(u::Float64)
    return 0.0 # Output
end

# In all these functions, ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
const fu = fu_pol           # Mapping function
const dfdu = dfdu_pol       # First derivative
const d2fdu2 = d2fdu2_pol   # Second derivative
const d3fdu3 = d3fdu3_pol   # Third derivative
const d4fdu4 = d4fdu4_pol   # Third derivative
