# Physical constants of the code
const G = 1.0 # Gravitational constant
const MTOTBATH = 1.0 # Total mass of the bath particles
# Caution : The total mass may be a little bit different
# and may vary along simulation (see Variables.jl)

# For GTE
const LAMBDA = 1.0 # Characteristic length-scale
const SIGMA = sqrt(LAMBDA*G*MTOTBATH) # Characteristic velocity-scale