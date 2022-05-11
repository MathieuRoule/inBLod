# Physical constants of the code
const G = 1.0 							# Gravitational constant
const Lambda = 1.0 						# Characteristic length-scale of the thermal equilibrium
const Mtot = 1.0 						# Total mass of the potential
const Sigma = sqrt(Mtot * G * Lambda)	# Characteristic velocity-scale of the thermal equilibrium
const mu = Mtot/NPART # Individual mass of the particles

# For other QSS with same total energy
if STATE == "Polytrope"
    const Alpha = 2.0 * Lambda
elseif STATE == "Plummer"
    const Alpha = 2.0 / pi * Lambda
end

# Julia constants
const NTHREADS = Threads.nthreads()

