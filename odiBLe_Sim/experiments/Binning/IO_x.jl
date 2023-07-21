##################################################
# At the end of all tries 
# We dump :
# - The positions histogram
##################################################
function dump_hist_x!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)                  # Opening the file with read and write
    #####
    write(file,"hist_x",hist_x)        # Dumping the table of position histogram
    #####
    close(file)
end


#################################################
# Saving the bins list/parameters at the end of the simulation
#################################################
function dump_bins_x_params!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)      # Opening the file with read and write
    #####
    # Positions bins
    #####
    write(file,"xlist",xbin_list)       # Dumping the bin centers list
    write(file,"nbins_x",NBINS_x)       # Dumping the nb of bins
    write(file,"dx",dx_bin)             # Dumping the bin width
    #####
    close(file)
end

#################################################
# Reading the bin list/parameters 
#################################################
const bins_x_params_list = ["xlist","nbins_x","dx"]
# Then use (from core/IO.jl) the function 
# att_read(namefile,bins_x_params_list)
