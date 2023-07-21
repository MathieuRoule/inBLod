##################################################
# At the end of each try 
# We dump :
# - The meanDE2 tab
# - The bincount tab (for further mean over try purpose)
##################################################
function dump_DE2_binned!(namefile::String,dataset::String,option::String="r+")
    #####
    file = h5open(namefile,option)                  # Opening the file with read and write
    #####
    write(file,dataset*"_meanDE2",meanDE2)          # Dumping the table of mean (\Delta E)^2
    write(file,dataset*"_bincount",bincount)        # Dumping the table of the number of particle in each bin
    write(file,dataset*"_meanEinbin",meanEinbin)    # Dumping the table of the mean energy in each bin
    #####
    close(file)
end


#################################################
# Saving the bins list/parameters at the end of the simulation
#################################################
function dump_bins_params!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)      # Opening the file with read and write
    #####
    write(file,"Elist",Ebin_list)       # Dumping the bin centers list
    write(file,"nbins",NBINS_E)
    write(file,"dE",dE_bin)             # Dumping the bin width
    write(file,"npart4DEE",NPART4DEE)   # Dumping the nb of particles used for statistics
    #####
    close(file)
end

#################################################
# Reading the bin list/parameters 
#################################################
const bins_params_list = ["Elist","nbins","dE","npart4DEE"]
# Then use (from core/IO.jl) the function 
# att_read(namefile,bins_params_list)