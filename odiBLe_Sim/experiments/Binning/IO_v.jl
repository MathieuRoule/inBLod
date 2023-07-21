##################################################
# At the end of all tries 
# We dump :
# - The velocities histogram
##################################################
function dump_hist_v!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)                  # Opening the file with read and write
    #####
    write(file,"hist_v",hist_v)        # Dumping the table of position histogram
    #####
    close(file)
end


#################################################
# Saving the bins list/parameters at the end of the simulation
#################################################
function dump_bins_v_params!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)      # Opening the file with read and write
    #####
    # Velocities bins
    #####
    write(file,"vlist",vbin_list)       # Dumping the bin centers list
    write(file,"nbins_v",NBINS_v)       # Dumping the nb of bins
    write(file,"dv",dv_bin)             # Dumping the bin width
    #####
    close(file)
end

#################################################
# Reading the bin list/parameters 
#################################################
const bins_v_params_list = ["vlist","nbins_v","dv"]
# Then use (from core/IO.jl) the function 
# att_read(namefile,bins_x_params_list)
