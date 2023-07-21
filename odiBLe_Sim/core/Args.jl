##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    #####
    # Number of particles
    #####
    "--nbathpart"
    help = "Number of bath particles"
    arg_type = Int64
    default = 1000
    "--ntestpart"
    help = "Number of test particles (saved at each test save)"
    arg_type = Int64
    default = 0
    #####
    # For integration
    #####
    "--scheme"
    help = "Integration scheme: symplectic/exact"
    arg_type = String
    default = "symplectic"
    "--order"
    help = "Order of the integration scheme (symplectic only): 2"
    arg_type = Int64
    default = 2
    "--dt"
    help = "Integration timestep (symplectic only)"
    arg_type = Float64
    default = 0.01
    #####
    # Warmup / simulation / saves times
    #####
    "--twarm"
    help = "Thermalisation time at the beginning before any save (not accounted for in the total time)"
    arg_type = Float64
    default = 0.0
    "--tmax"
    help = "Total integration time"
    arg_type = Float64
    default = 100.0
    "--tsavescale"
    help = "Time scale for saves (linear/log2/log10)"
    arg_type = String
    default = "linear"
    "--tensave"
    help = "Save rate in time between 2 saves (only total energy and momentum)"
    arg_type = Float64
    default = 1.0
    "--ttestsave"
    help = "Save rate in time bewteen 2 saves of test particles (tabx and tabv)"
    arg_type = Float64
    default = 50.0
    "--tfullsave"
    help = "Save rate in time bewteen 2 full saves (tabx and tabv)"
    arg_type = Float64
    default = 100.0
    #####
    # Random seed (for reproductibility purpose)
    #####
    "--seed"
    help = "Random seed"
    arg_type = Int64
    default = 1
    #####
    # Initial sampling
    #####
    "--fbathsampling"
    help = "Initial positions/velocities sampling for bath particles : GTE/Plummer/Polytrope/rand/bar"
    arg_type = String
    default = "GTE"
    "--ftestsampling"
    help = "Initial positions/velocities sampling for test particles : GTE/Plummer/Polytrope/rand/unif/bar/circle/oneplusGTE"
    arg_type = String
    default = "GTE"
    "--msampling"
    help = "Initial mass sampling : equal/testnomassp/bigp/bigandnomass"
    arg_type = String
    default = "equal"
    #####
    # Landau or Balescu-Lenard type experiment (Landau --> bath particles follow the mean potential)
    #####
    "--coupling"
    help = "Balescu-Lenard (All part. follow the inst. pot.) or Landau (bath part. follow the mean pot.) or Forcing_BL (massive 'test' part. follow the mean 
    pot. and bath part. the inst. pot.) or Forcing_Landau (first massive test part. and bath follow the mean pot. and other massless test part. the inst. pot.)"
    arg_type = String
    default = "Balescu-Lenard"
    #####
    # Number of different random initial sampling to do
    #####
    "--ntry"
    help = "Number of try (different initial sampling)"
    arg_type = Int64
    default = 1
    #####
    # Outputs
    #####
    "--filename"
    help = "Path + filename to use to dump infos into a HDF5 file (without .hf5)"
    arg_type = String
    default = "./dump"
end
parsed_args = parse_args(tabargs)
