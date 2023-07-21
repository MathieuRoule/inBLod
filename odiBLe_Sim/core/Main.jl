include("Packages.jl") # To load all the needed packages
##################################################
const to = TimerOutput() # Create the timer object
##################################################
include("Args.jl") # Parsing the command-line
##################################################
# General parameters
##################################################

const NBATHPART = parsed_args["nbathpart"] # Number of bath particles
const NTESTPART = parsed_args["ntestpart"] # Number of test particles 
const NPART = NBATHPART + NTESTPART # Total number of particles (bath + test)

#####
const SCHEME = parsed_args["scheme"] # Determining the integration scheme
const ORDER = parsed_args["order"] # Order of the integration scheme. ATTENTION, symplectic only

const DT = parsed_args["dt"] # Integration timestep. ATTENTION, symplectic only.
const TWARM = parsed_args["twarm"] # Thermalisation/warmup time (no save during this time)
const TMAX = parsed_args["tmax"] # Total integration time (without warmup time)
const TSAVESCALE = parsed_args["tsavescale"] # Time scale to use for saves (linear, 2-logarithmic)
const TENSAVE = parsed_args["tensave"] # Save rate in time bewteen 2 saves (only total energy and momentum)
const TTESTSAVE = parsed_args["ttestsave"] # Save rate in time bewteen 2  part saves (tabx and tabv of few particles)
const TFULLSAVE = parsed_args["tfullsave"] # Save rate in time bewteen 2 full saves (tabx and tabv)

const SEED = parsed_args["seed"] # Random seed

const INITFBATH = parsed_args["fbathsampling"] # Initial positions/velocities sampling of bath particles
const INITFTEST = parsed_args["ftestsampling"] # Initial positions/velocities sampling of test particles
const INITM = parsed_args["msampling"] # Initial mass sampling (only 1 function for both bath and test particles)
const COUPLING = parsed_args["coupling"] # "Balescu-Lenard" or "Landau"

#####
const NTRY = parsed_args["ntry"]
const FILENAME = parsed_args["filename"]


##################################################
# Sanity checks on the command-line arguments
##################################################

#####
# Check on the scheme
#####
if ((SCHEME != "symplectic") && (SCHEME != "exact"))
    error("ERROR: UNKNOWN SCHEME") # Unknown integration scheme
end
#####
# Check on the order
#####
if (ORDER != 2)
    error("ERROR: UNKNOWN ORDER") # Unknown integratin order
end

#####
# Check on the time save scale
#####
if ( (TSAVESCALE != "linear") && (TSAVESCALE != "log2") && (TSAVESCALE != "log10") )
    error("ERROR: UNKNOWN TIME SAVE SCALE") # Unknown time save scale
end

#####
# Check on the position/velocities initialization for both bath and test particles
#####
if ((INITFBATH != "GTE") && (INITFBATH != "Polytrope") && (INITFBATH != "Plummer") && (INITFBATH != "rand") && (INITFBATH != "bar"))
    error("ERROR: UNKNOWN BATH INITIALIZATION SAMPLING") # Unknown sampling
end
if ((INITFTEST != "GTE") && (INITFBATH != "Polytrope") && (INITFBATH != "Plummer") && (INITFTEST != "rand") && (INITFTEST != "unif") && (INITFTEST != "bar") && (INITFTEST != "circle") && (INITFTEST != "oneplusGTE"))
    error("ERROR: UNKNOWN TEST INITIALIZATION SAMPLING") # Unknown sampling
end
#####
# Check on the mass initialization
#####
if ((INITM != "equal") && (INITM != "testnomassp") && (INITM != "bigp") && (INITM != "bigandnomass"))
    error("ERROR: UNKNOWN INITIALIZATION MASS SAMPLING") # Unknown sampling
end

#####
# Check on the coupling
#####
if ((COUPLING != "Balescu-Lenard") && (COUPLING != "Landau") && (COUPLING != "Forcing_BL") && (COUPLING != "Forcing_Landau"))
    error("ERROR: UNKNOWN COUPLING") # Unknown coupling
end

##################################################

include("Constants.jl") # Physical constants
include("Variables.jl") # Declaration of variables and arrays
include("IO.jl") # To be able to dump/read .hf5 files
include("Sampling.jl") # To sample the initial conditions
include("Energy.jl") # For energy/momentum calculation
include("Forces.jl") # Forces computation (Landau or Balescu-Lenard)



####
if (SCHEME == "symplectic") # We load the symplectic integration functions
    include("Symplectic.jl") # Symplectic integration
end
#####
if (SCHEME == "exact") # We load the exact integration functions
    include("Exact.jl") # Exact integration
end
#####
