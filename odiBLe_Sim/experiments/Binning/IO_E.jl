##################################################
# At the end of all tries 
# We dump :
# - The energies histogram
##################################################
function dump_hist_E!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)                  # Opening the file with read and write
    #####
    write(file,"hist_E",hist_E)        # Dumping the table of energy histogram
    #####
    close(file)
end


#################################################
# Saving the bins list/parameters at the end of the simulation
#################################################
function dump_bins_E_params!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)      # Opening the file with read and write
    #####
    # Energy bins
    #####
    write(file,"Elist",Ebin_list)       # Dumping the bin centers list
    write(file,"nbins_E",NBINS_E)       # Dumping the nb of bins
    write(file,"dE",dE_bin)             # Dumping the bin width
    #####
    close(file)
end

#################################################
# Reading the bin list/parameters 
#################################################
const bins_E_params_list = ["Elist","nbins_E","dE"]
# Then use (from core/IO.jl) the function 
# att_read(namefile,bins_E_params_list)
