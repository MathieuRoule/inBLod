include("./../../core/Main.jl") # Loading all the code

include("./Savedqt_binDEE.jl")    # Creating arrays for saves
include("./DEE_IO.jl")            # Dumping function 
include("./DEE_Binning.jl")       # (\Delta E)^2 mean calculation on each bin

###################################################
# Print advertisers
###################################################
# For Landau like experiment
if COUPLING == "Landau"
    println("WARNING : Landau mode of forces computation activated")
end


###################################################
# Saving function along simulation needed 
# for ../core/Run.jl Run!() function.
###################################################
function saves_at_en_rate!(ntry::Int64,count::Int64)
    save_en!(count)
end

function saves_at_test_rate!(ntry::Int64,count::Int64)
    DE2_calc!(count)       # Calculation of the mean (\Delta E)^2
    save_current_time!(count,testsaved_times)   # Saving time
end

function saves_at_full_rate!(ntry::Int64,count::Int64)
end

function init_saves!(ntry::Int64)
    #####
    # Initialization of the bin
    #####
    init_binning!()
    #####
    saves_at_en_rate!(ntry,1)   # Initial energy save
    saves_at_test_rate!(ntry,1) # Initial histogram save
end

function saves_at_endtry!(ntry::Int64)
    dump_endtry!(namefile,"test_"*string(ntry))
    dump_DE2_binned!(namefile,"test_"*string(ntry))
end

function saves_at_end!()
    att_dump!(namefile)             # Simulation general parameters
    dump_bins_params!(namefile)     # DE2 specific parameters

end

###################################################
# Run
###################################################
include("./../../core/Run.jl")
Run!()