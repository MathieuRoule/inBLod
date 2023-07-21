################################################
# Generic Julia script to run a simulation with 
# the given parameters from core/Args.jl.
#
# Informations to save from the simulation 
# at 3 differents rates (en, test, full)
# have to be defined according to the experiment (see Run! function infos)
# and then this program can be called to perform the simulations.
#
# CAUTION : "core/Main.jl" has to be included before
################################################


########################
# File creation
########################
if !(@isdefined namefile)
    error("namefile is not defined: core/Main.jl has to be included before this file")
end
file = h5open(namefile,"w")                     # File creation, hence the option w.
close(file)

########################
# Fill all the array with zeros (at the beginning of a new try)
# To be sure there is no influence of a try over the following one.
# See all the array in "./Variables.jl"
########################
function fillingzeros!()
    TIME[1] = 0.0
    MTOT[1] = 0.0
    #####
    v_mf[1] = 0.0
    #####
    fill!(tabm,0.0)
    fill!(tabx,0.0)
    fill!(tabv,0.0)
    fill!(tabvdot,0.0)
    fill!(tabperm,0.0)
    for i = 1:NPART
        tabperm[i] = i # Initialization of the permutation array so it is a permutation of 1:N from the beginning.
    end
    fill!(tabcumm,0.0)
    fill!(tabcummx,0.0)
    #####
    fill!(ensaved_times,0.0)
    fill!(tabE,0.0)
    fill!(tabP,0.0)
    #####
    fill!(testsaved_times,0.0)
    #####
    fill!(fullsaved_times,0.0)
end


########################
# Initialization and warm up (during which no save is done).
########################
function initialization_warmup!()
    fillingzeros!() # Re-initialization of the parameters at the beginning of the try
    #####
    # Initial sampling (before warm up)
    #####
    init_bathsampling!()   # Performing the initial sampling of the bath particles
    init_testsampling_1!() # Performing the initial sampling of the test particles
    #####
    # Defining the mass of particles (for the warm up)
    #####
    init_mass_1!() 
    #####
    # Centering distribution (frame changing so the total momentum is 0)
    #####
    centered_frame!()

    #####
    # Warm-up (time not accounted for in the total time)
    #####
    while TIME[1]<TWARM
        @timeit to "warmupintegration" integrate_DT!()
    end
    TIME[1] = 0.0
    #####
    # After warm up initializations ((m,x,v) of test particles)
    #####
    init_testsampling_2!()
    init_mass_2!() 
    #####
    # Centering distribution (frame changing so the total momentum is 0)
    #####
    centered_frame!()

end

########################
# Run function
#
# To specify before use :
#
# - init_saves!(ntry) function -> to do right after the initializations and warm up
#
# - saves_at_en_rate!(ntry,count) function
# - saves_at_test_rate!(ntry,count) function
# - saves_at_full_rate!(ntry,count) function
#
# - saves_at endtry!(ntry) function
# - saves_at_end!() function
########################
function Run!()
    
    ##################################################
    # Integration for NSTEPS timesteps
    ##################################################
    Random.seed!(SEED) # Fixing the random seed
    ##################################################

    for nt = 1:NTRY
        @timeit to "print" println("Try number ",nt," over ",NTRY)
        
        #####
        # (m,x,v) initialization and warmup for test and bath particles
        #####
        @timeit to "init" initialization_warmup!()
        
        #####
        # Initial saves
        #####
        @timeit to "save" init_saves!(nt)
        
        #####
        # Counters
        #####
        ensave_count = 1
        testsave_count = 1
        fullsave_count = 1
        tes = TENSAVE
        tts = TTESTSAVE
        tfs = TFULLSAVE

        #####
        # Integration
        #####
        while TIME[1]<TMAX
            @timeit to "integration" integrate_DT!()
            #####
            # Energy save
            #####
            if TIME[1] >= tes
                ensave_count += 1
                @timeit to "save" saves_at_en_rate!(nt,ensave_count) 
                if TSAVESCALE == "linear"
                    tes += TENSAVE
                elseif TSAVESCALE == "log2"
                    tes *= 2.0
                elseif TSAVESCALE == "log10"
                    tes *= 10.0
                end
            end
            #####
            # Part save
            #####
            if (TIME[1] >= tts)
                testsave_count += 1
                @timeit to "save" saves_at_test_rate!(nt,testsave_count)
                if TSAVESCALE == "linear"
                    tts += TTESTSAVE
                elseif TSAVESCALE == "log2"
                    tts *= 2.0
                elseif TSAVESCALE == "log10"
                    tts *= 10.0
                end
            end
            #####
            # Full save
            #####
            if TIME[1] >= tfs
                fullsave_count += 1
                @timeit to "save" saves_at_full_rate!(nt,fullsave_count)
                if TSAVESCALE == "linear"
                    tfs += TFULLSAVE
                elseif TSAVESCALE == "log2"
                    tfs *= 2.0
                elseif TSAVESCALE == "log10"
                    tfs *= 10.0
                end
            end
        end
        @timeit to "save" saves_at_endtry!(nt)
    end
    @timeit to "save" saves_at_end!()
        
    show(to)
    println("\n END")
end