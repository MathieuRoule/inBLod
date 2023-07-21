###########################
# This program return the contribution 
# of each resonance (n,np) to the 
# flux, diffusion and friction 
# at a given ra.
###########################

# Include the main computation code
include("./../../core/Main.jl")

# Packages needed for this experiment
pkgcheck("TimerOutputs","ProgressMeter")

using ProgressMeter
using TimerOutputs
const to = TimerOutput() # Create the timer object


# Creation and initialization of the structures needed for computation
@timeit to "init" include("./../../core/Structscreateinit.jl")

#####
# Considered ra
#####
if (STATE=="GTE")
    const ra = Lambda
else
    const ra = Alpha
end
const val_Omega = Omegafromra(ra) # Current frequency

const tab_AJ_contrib    = zeros(Float64,NMAX,NMAX)
const tab_DJJ_contrib   = zeros(Float64,NMAX,NMAX)
const tab_Flux_contrib  = zeros(Float64,NMAX,NMAX)

######
function rescontrib!()

    p = Progress(NMAX^2,desc="Computations: ")

    Threads.@threads for n=1:NMAX
	k = Threads.threadid() # Thread in charge
        ndotomega = n*val_Omega # Current resonance frequency
        for np=1:NMAX # Loop over np. ATTENTION, the `step' is 2
            rap = get_Resonance(ndotomega,np) # Solving for the resonance condition
            #####
            if (rap != nothing) # We have found a solution to the resonance condition
                (tab_AJ_contrib[n,np], tab_DJJ_contrib[n,np], tab_Flux_contrib[n,np]) = get_resonance_contribution(n,np,ra,rap,ndotomega,structMultipole[k],structBasis[k],structBL[k]) # Contribution from the current resonance
            end
            #####
            next!(p)
        end
    end
end
#####
@timeit to "calcul" rescontrib!()

function dump_rescontrib!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option) # Opening the file with write-only, i.e. it gets overwritten
    #####
    # Simulation results
    #####
    write(file,"considered_ra",ra)
    write(file,"tabAJ",tab_AJ_contrib) 
    write(file,"tabDJJ",tab_DJJ_contrib)  
    write(file,"tabFlux",tab_Flux_contrib) 
    #####
    close(file)
end

att_dump!(filename,"w")
dump_rescontrib!(filename)

show(to)
println("\n END")
