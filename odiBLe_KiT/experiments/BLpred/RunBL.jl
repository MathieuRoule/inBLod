###########################
# Return the total
# flux, diffusion and friction 
# at different ra in a HDF5 file.
###########################

try using TimerOutputs catch; println("The TimerOutputs package is needed."); Pkg.add("TimerOutputs"); using TimerOutputs end
const to = TimerOutput() # Create the timer object

# Include the main computation code
include("./../../core/Main.jl")
# Creation and initialization of the structures needed for computation
@timeit to "init" include("./../../core/Structscreateinit.jl")


const tab_ra    = collect(0.01:0.05:5.01)
const nb_ra     = length(tab_ra)
const tab_E     = zeros(Float64,nb_ra)
const tab_J     = zeros(Float64,nb_ra)
const tab_Om    = zeros(Float64,nb_ra)
const tab_domdr = zeros(Float64,nb_ra)
const tab_dJdr  = zeros(Float64,nb_ra)
const tab_F0J   = zeros(Float64,nb_ra)
const tab_dF0dJ = zeros(Float64,nb_ra)

const tab_AJ    = zeros(Float64,nb_ra)
const tab_DJJ   = zeros(Float64,nb_ra)
const tab_DEE   = zeros(Float64,nb_ra)
const tab_Flux 	= zeros(Float64,nb_ra)

######
function tabdiffcoef!()
    Threads.@threads for i=1:nb_ra
        k               = Threads.threadid()
        ra              = tab_ra[i]
        tab_E[i]        = psi(ra)
        tab_J[i]        = Jfromra(ra)
        tab_Om[i]       = Omegafromra(ra)
        tab_domdr[i]    = dOmegadra(ra)
        tab_dJdr[i]     = JacJtora(ra)
        tab_F0J[i]      = DF(ra)
        tab_dF0dJ[i]    = dDFdJ(ra)
        (tab_AJ[i], tab_DJJ[i], tab_Flux[i]) = get_AJ_DJJ_Flux(ra,structMultipole[k],structBasis[k],structBL[k])
        tab_DEE[i]      = get_DEE_fromOm_DJJ(tab_Om[i],tab_DJJ[i])
        println(i)
    end
end
#####
@timeit to "calcul" tabdiffcoef!()

function dump_BLpred!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option) # Opening the file with write-only, i.e. it gets overwritten
    #####
    # Simulation results
    #####
    write(file,"tabra",tab_ra) 
    write(file,"tabE",tab_E) 
    write(file,"tabJ",tab_J) 
    write(file,"tabOm",tab_Om) 
    write(file,"tabdomdr",tab_domdr) 
    write(file,"tabdJdra",tab_dJdr) 
    write(file,"tabF0J",tab_F0J) 
    write(file,"tabdF0dJ",tab_dF0dJ)

    write(file,"tabAJ",tab_AJ) 
    write(file,"tabDJJ",tab_DJJ) 
    write(file,"tabDEE",tab_DEE) 
    write(file,"tabFlux",tab_Flux) 
    #####
    close(file)
end

att_dump!(filename,"w")
dump_BLpred!(filename)

show(to)
println("\n END")
