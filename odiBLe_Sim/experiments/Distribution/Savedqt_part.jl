##################################################
# Quantities calculated/saved at TENSAVE rate
##################################################

# Nothing more than Etot and Ptot (already in core)

##################################################
# Quantities calculated/saved at TTESTSAVE rate
# We recall the number of test particles saves NTESTSAVES (in core)
##################################################
@timeit to "allocation" const tabtestsavex = zeros(Float64,NTESTPART,NTESTSAVES)   # Array used to store the part save (only one dataset in the HDF5 file)
@timeit to "allocation" const tabtestsavev = zeros(Float64,NTESTPART,NTESTSAVES)   # Array used to store the part save (only one dataset in the HDF5 file)
##################################################
# Quantities calculated/saved at TFULLSAVE rate (except full state - tabx and tabv - already in core)
# We recall the number of full saves is NFULLSAVES (in core)
##################################################

# Nothing more than tabx and tabv (already in core)