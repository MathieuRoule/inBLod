##################################################
# Sinusoidal Effective anomaly
##################################################
function fu_sin(u::Float64)
    return sin(0.5 * pi * u) # Output
end
##########
# Derivative of the sinusoidal effective anomaly
##########
function dfdu_sin(u::Float64)
    return 0.5 * pi * cos(0.5 * pi * u)  # Output
end
##########
# Second derivative of the sinusoidal effective anomaly
##########
function d2fdu2_sin(u::Float64)
    return -(0.5 * pi)^(2) * sin(0.5 * pi * u)  # Output
end
##########
# Third derivative of the sinusoidal effective anomaly
##########
function d3fdu3_sin(u::Float64)
    return -(0.5 * pi)^(3) * cos(0.5 * pi * u)  # Output
end
##########
# Fourth derivative of the sinusoidal effective anomaly
##########
function d4fdu4_sin(u::Float64)
    return (0.5 * pi)^(4) * sin(0.5 * pi * u)  # Output
end

# In all these functions, ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
const fu = fu_sin           # Mapping function
const dfdu = dfdu_sin       # First derivative
const d2fdu2 = d2fdu2_sin   # Second derivative
const d3fdu3 = d3fdu3_sin   # Third derivative
const d4fdu4 = d4fdu4_sin   # Third derivative