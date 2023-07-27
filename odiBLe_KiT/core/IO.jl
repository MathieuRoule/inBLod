###############################################
# Attributes dump at the end of the program execution
# Datasets containing attributes caracterizing the computation
###############################################
function att_dump!(namefile::String,option::String="r+")
    file = h5open(namefile,option)
    #####
    write(file,"npart",NPART)
    write(file,"nmax",NMAX)
    #####
    write(file,"state",STATE)
    #####
    write(file,"mapping",MAPPING)
    write(file,"K_MF",K_MF)
    #####
    write(file,"coupling",COUPLING)
    write(file,"K_FT",K_FT)
    #####
    write(file,"basis",BASIS)
    write(file,"nbasiselt",NBASISELMT)
    write(file,"L_period",L_p)
    #####
    write(file,"BLmethod",BLMETHOD)
    write(file,"nmax_resmat",NMAXRESMAT)
    write(file,"BL_ra_max",BL_RA_MAX)
    write(file,"K_BL",K_BL)
    #####
    close(file)
end
