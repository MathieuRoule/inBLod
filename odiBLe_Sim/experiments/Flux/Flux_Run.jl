include("./../../core/Main.jl") # Loading all the code

include("./../Binning/hist_E.jl")     # Histogram in energy handling

###################################################
# Print advertisers
###################################################
# For experiments which do not correspond to wake calculation
if ((COUPLING != "Balescu-Lenard") && (COUPLING != "Landau"))
    println("WARNING : method for forces computation different from 'Balescu-Lenard' or 'Landau' activated")
end

###################################################
# Saving function along simulation needed 
# for ../core/Run.jl Run!() function.
###################################################
function saves_at_en_rate!(ntry::Int64,count::Int64)
    save_en!(count)
end

function saves_at_test_rate!(ntry::Int64,count::Int64) 
    hist_E_calc!(count)                         # Histograms save
    save_current_time!(count,testsaved_times)   # Saving time
end

function saves_at_full_rate!(ntry::Int64,count::Int64)
end

function init_saves!(ntry::Int64)
    saves_at_en_rate!(ntry,1)   # Initial energy save
    saves_at_test_rate!(ntry,1) # Initial histogram save
end

function saves_at_endtry!(ntry::Int64)
    dump_endtry!(namefile,"test_"*string(ntry))
end

function saves_at_end!()
    dump_hist_E!(namefile)          # Dumping histograms
    dump_bins_E_params!(namefile)   # Dumping histograms parameters
    att_dump!(namefile)             # Dumping all the simulation parameters
end

###################################################
# Run
###################################################
include("./../../core/Run.jl")
Run!()