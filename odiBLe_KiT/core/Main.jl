##################################################
include("Packages.jl") 		# Loading the needed packages
include("Args.jl") 			# Parsing the command-line
include("Constants.jl") 	# Physical constants
include("Mean.jl") 			# Description of the mean system
include("Resonances.jl") 	# To find the resonant locations
include("PsiSQ.jl") 		# To compute the coupling coefficients
include("Flux.jl") 			# To compute the friction, diffusion coefficients and the total flux
include("IO.jl")			# Dump the results
##################################################
