########################
# Dumping file path
########################
if (@isdefined FILENAME)
    # Filename for the HDF5 dump (must contain the relative path from where the program was launched)
    const namefile = FILENAME*".hf5"
end
##################################################
# Function to dump the current state
# of the simulation
# Inputs:
# + namefile: name of the file
# + dataset: prefix name for the dataset
# + option for the reading: by default, it is not overwritten
##################################################
function dump_fullstate!(namefile::String,dataset::String,k::Int64,option::String="r+")
    #####
    file = h5open(namefile,option) # Opening the file with write-only, i.e. it gets overwritten
    #####
    write(file,dataset*"_fullsaved_x",tabx) # Dumping the table of positions (mean DF reference frame = simulation RF)
    write(file,dataset*"_fullsaved_v",map(v_from_sim_to_meanDF_RF,tabv)) # Dumping the table of velocities in the mean DF reference frame
    fullsaved_times[k]=TIME[1]
    #####
    close(file)
end

##################################################
# Saving current time in a time save array
# (ensaved_times, testsaved_times or fullsaved_times)
# 2 arguments expected :
#   - n -> the position in the array (i.e. number of saves counter)
#   - timetab -> which 
##################################################
function save_current_time!(n::Int64,timetab::Array{Float64})
    timetab[n] = TIME[1]
end

##################################################
# At the end of each try 
# We dump :
# - The momentum center of the realization
# - The total energy and the momentum along simulation
# - The times when energy/momentum, test and full save were made
##################################################
function dump_endtry!(namefile::String,dataset::String,option::String="r+")
    #####
    file = h5open(namefile,option) # Opening the file with read and write
    #####
    write(file,dataset*"_vmf",v_mf[1]) # Dumping the momentum recentering of the realization
    #####
    write(file,dataset*"_tabE",tabE) # Dumping the table of total energy along the simulation
    write(file,dataset*"_tabP",tabP) # Dumping the table of total momentum along the simulation
    #####
    write(file,dataset*"_ensaved_t",ensaved_times)     # Dumping the table of times where energy and momentum where saved
    write(file,dataset*"_testsaved_t",testsaved_times) # Dumping the table of times where part state where saved (tabx and tabv)
    write(file,dataset*"_fullsaved_t",fullsaved_times) # Dumping the table of times where full state where saved (tabx and tabv)
    #####
    close(file)
end

###############################################
# Attributes dump after all try (at the end of the program execution)
# Datasets containing attributes caracterizing the simulation
###############################################
function att_dump!(namefile::String,option::String="r+")
    file = h5open(namefile,option)
    #####
    write(file,"nbathpart",NBATHPART)
    write(file,"ntestpart",NTESTPART)
    #####
    write(file,"npart",NPART)
    #####
    write(file,"scheme",SCHEME)
    write(file,"order",ORDER)
    write(file,"dt",DT)
    #####
    write(file,"twarm",TWARM)
    write(file,"tmax",TMAX)
    write(file,"tsavescale",TSAVESCALE)
    write(file,"tensave",TENSAVE)
    write(file,"ttestsave",TTESTSAVE)
    write(file,"tfullsave",TFULLSAVE)
    #####
    write(file,"nensaves",NENSAVES)
    write(file,"ntestsaves",NTESTSAVES)
    write(file,"nfullsaves",NFULLSAVES)
    #####
    write(file,"seed",SEED)
    #####
    write(file,"fbathsampling",INITFBATH)
    write(file,"ftestsampling",INITFTEST)
    write(file,"msampling",INITM)
    #####
    write(file,"coupling",COUPLING)
    #####
    write(file,"ntry",NTRY)
    #####
    close(file)
end

const att_list = ["nbathpart","ntestpart","npart","scheme","order","dt","twarm","tmax","tsavescale",
                  "tensave","ttestsave","tfullsave","nensaves","ntestsaves","nfullsaves",
                  "seed","fbathsampling","ftestsampling","msampling","coupling","ntry"]

###############################################
# Attributes read
# Return a Dict object containing all the 
# attributes values of the file
###############################################
function att_read(namefile::String,attributes::Array{String,1}=att_list)
    att_dict = Dict{String,Any}()   # Attributes dictionnary
    #####
    file = h5open(namefile,"r")     # read only
    #####
    for attribute in attributes
        setindex!(att_dict,read(file,attribute),attribute)
    end
    #####
    close(file)
    #####
    return att_dict
end

###############################################
# Attributes dump from the Dict object returned
# by att_read() function
###############################################
function att_dump_fromdict!(namefile::String,option::String="r+",dict::Dict{String,Any}=att_list)
    file = h5open(namefile,option)
    #####
    for (key, value) in dict
        write(file,key,value)
    end
    #####
    close(file)
end