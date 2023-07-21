include("./../../core/Main.jl") # Loading all the core code

include("Savedqt_part.jl") # Creating arrays for test saves
include("IO_part.jl")      # Specific dumping functions

###################################################
# Print advertisers
###################################################
# For Landau like experiment
if COUPLING == "Landau"
    @timeit to "print" println("WARNING : Landau mode of forces computation activated")
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
    # Initial fullsave
    # Dumping the initial conditions. 
    #####
    dump_fullstate!(namefile,"test_"*string(ntry)*"_t_1",1)
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
    dump_fullstate!(namefile,"test_"*string(ntry)*"_t_"*string(count),count)
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