##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    #####
    #   General Part
    #####
    "--npart"
    help = "Number of particles"
    arg_type = Int64
    default = 10000
    "--nmax"
    help = "Maximum resonance number"
    arg_type = Int64
    default = 50
    #####
    #   Considered Quasi-Stationary state 
    #####
    "--state"
    help = "QSS considered: GTE/Plummer/Polytrope"
    arg_type = String
    default = "GTE"
    #####
    #   For mean functions calcultation (J(ra), Omega(ra)...)
    #####
    "--mapping"
    help = "Mapping used: polynomial/sinusoidal"
    arg_type = String
    default = "polynomial"
    "--K_MF"
    help = "Number of intervals in u (effective anomaly) for mean function integration (on [0,1])"
    arg_type = Int64
    default = 100
    #####
    #   Bare or dressed susceptibility coefficient / Multipole (only for bare) or Basis choice
    #####
    "--coupling"
    help = "Coupling considered: Balescu-Lenard/BL-PdS/Landau-Basis/Landau-Multipole"
    arg_type = String
    default = "Landau-Basis"
    "--K_FT"
    help = "Number of nodes in u (effective anomaly) for FT of basis elements (or Multipole expansion)"
    arg_type = Int64
    default = 1000
    #####
    #   If Basis :
    #####
    "--basis"
    help = "Considered basis: only CosSin available at the moment"
    arg_type = String
    default = "CosSin"
    "--nbasiselt"
    help = "Number of basis elements (*2 - odd and even)"
    arg_type = Int64
    default = 50
    "--L_period"
    help = "Characteristic length of periodisation of the potential (*Lambda)"
    arg_type = Float64
    default = 10.0
    #####
    #   If dressed, choices for Response Matrix calculation
    #####
    "--BLmethod"
    help = "If Balescu-Lenard coupling, matrix method considered: Splitting/Legendre"
    arg_type = String
    default = "Legendre"
    "--nmax_resmat"
    help = "If Balescu-Lenard coupling, maximum resonance number for the response matrix calculation"
    arg_type = Int64
    default = 20
    "--BL_ra_max"
    help = "ra_max for integration in response matrix calculation (*Lambda)"
    arg_type = Float64
    default = 5.0
    "--K_BL"
    help = "Number of splits in splitting method or interpolation nodes in Legendre method (response matrix)"
    arg_type = Int64
    default = 100
    #####
    #   For outputs
    #####
    "--filename"
    help = "Path + name of the HDF5 file that will be created"
    arg_type = String
    default = "./KinTh"
end
parsed_args = parse_args(tabargs)

##################################################
const NPART    = parsed_args["npart"]               # Number of particles
const NMAX     = parsed_args["nmax"]                # Maximum resonance number
##########
const STATE    = parsed_args["state"]               # Considered quasi-stationary state
##########
const MAPPING  = parsed_args["mapping"]             # Mapping function used

const K_MF   = parsed_args["K_MF"]                  # Nodes in u for smooth function integration
##########
const COUPLING = parsed_args["coupling"]            # Coupling considered
const K_FT      	= parsed_args["K_FT"]           # Nodes in u (effective anomaly) for FT
##########
const BASIS         = parsed_args["basis"]          # Basis elements to use
const NBASISELMT    = parsed_args["nbasiselt"]      # Number of basis elts to use (in fact NBASISELMT odd elts +  NBASISELMT even elts --> 2*NBASISELMT)
const L_p           = parsed_args["L_period"]       # For SinCos Basis, length of periodization of the potential (to multiply by Lambda)
##########
const BLMETHOD      = parsed_args["BLmethod"]       # Response matrix method considered
const NMAXRESMAT    = parsed_args["nmax_resmat"]    # Maximum resonance number for response matrix
const BL_RA_MAX     = parsed_args["BL_ra_max"]      # Maximal ra of integration for response matrix (to multiply by Lambda)
const K_BL          = parsed_args["K_BL"]           # Number of splits in splitting method or interpolation nodes in Legendre method 
##########
const filename = parsed_args["filename"]*".hf5"     # HDF5 file name


##################################################
# Sanity checks on the command-line arguments
##################################################

#####
# Check on the mappings
#####
if ((MAPPING != "polynomial") && (MAPPING != "sinusoidal"))
    error("ERROR: UNKNOWN MAPPING") # Unknown choice for the mapping function
end

#####
# Check on the QSS
#####
if ((STATE != "GTE") && (STATE != "Plummer") && (STATE != "Polytrope"))
    error("ERROR: UNKNOWN STATE") # Unknown choice for the quasi-stationary state
end

#####
# Check on the couplings
#####
# Balescu-Lenard: Basis method + Response matrix
# BL-PdS ("Pain de Sucre"): Balescu-Lenard completed by undressed at low separation (using Landau-Multipole - Landau-Basis)
# Landau-Basis: Basis method + NO collective effects
# Landau-Multipole: Multipole method + NO collective effects
if ((COUPLING != "Balescu-Lenard") && (COUPLING != "BL-PdS") && (COUPLING != "Landau-Basis") && (COUPLING != "Landau-Multipole"))
    error("ERROR: UNKNOWN COUPLING") # Unknown choice for the coupling
end

#####
# Check on the basis
#####
if ((BASIS != "CosSin"))
    error("ERROR: UNKNOWN BASIS") # Unknown choice for the basis functions
end

#####
# Check on the response matrix method
#####
if ((BLMETHOD != "Splitting") && (BLMETHOD != "Legendre"))
    error("ERROR: UNKNOWN BLMETHOD") # Unknown choice for the Response matrix method
end