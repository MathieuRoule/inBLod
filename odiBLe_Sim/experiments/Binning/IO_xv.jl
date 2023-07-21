##################################################
# At the end of all tries 
# We dump :
# - The positions/velocities histogram
##################################################
function dump_hist_xv!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)                  # Opening the file with read and write
    #####
    write(file,"hist_xv",hist_xv)      # Dumping the table of position/velocity histogram
    #####
    close(file)
end


#################################################
# Saving the bins list/parameters at the end of the simulation
#################################################
function dump_bins_xv_params!(namefile::String,option::String="r+")
    #####
    file = h5open(namefile,option)      # Opening the file with read and write
    #####
    # Position bins
    #####
    write(file,"xlist",xbin_list)       # Dumping the bin centers list
    write(file,"nbins_x",NBINS_x)       # Dumping the nb of bins
    write(file,"dx",dx_bin)             # Dumping the bin width
    #####
    # Velocity bins
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
const bins_xv_params_list = ["xlist","nbins_x","dx","vlist","nbins_v","dv"]
# Then use (from core/IO.jl) the function 
# att_read(namefile,bins_xv_params_list)
