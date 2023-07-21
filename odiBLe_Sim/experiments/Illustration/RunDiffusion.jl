include("./../../core/Main.jl") # Loading all the code

include("./../Distribution/Savedqt_part.jl") # Creating arrays for test saves
include("./../Distribution/IO_part.jl")      # Specific dumping functions

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
function init_saves!(ntry::Int64)
    #####
    # Initial calculation of the energy
    #####
    save_en!(1)
    #####
    # Initial part save
    #####
    if (NTESTPART != 0)
        save_teststate!(1) 
    end
end

function saves_at_en_rate!(ntry::Int64,count::Int64)
    save_en!(count)
end

function saves_at_test_rate!(ntry::Int64,count::Int64)
    if NTESTPART != 0 
        save_teststate!(count)
    end
end

function saves_at_full_rate!(ntry::Int64,count::Int64)
end

function saves_at_endtry!(ntry::Int64)
    dump_endtry!(namefile,"test_"*string(ntry))
    dump_endtry_testpart!(namefile,"test_"*string(ntry))
end

function saves_at_end!()
    att_dump!(namefile)
end

###################################################
# Run
###################################################
include("./../../core/Run.jl")
Run!()