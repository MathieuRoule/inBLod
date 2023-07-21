##################################################
# This generic program merge the results dumped 
# in differents files into 1 sumup file.
#
# The files' names have to follow the syntax
# "../data/*i.jl" with i in range(imin,stop=imax).
#
# Specific merging function depending on the 
# content dumped can be created for each experience
# and used through the Sumup function (see below).
##################################################

using ArgParse # To be able to parse command-line arguments
using HDF5 # To be able to dump/read HDF5 files

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--filenames"
    help = "Relative path + Beginning of the names of the files to read (used for sumup file name)"
    arg_type = String
    default = "dump"
    "--fileidmin"
    help = "Minimum file id to read"
    arg_type = Int64
    default = 1
    "--fileidmax"
    help = "Maximum file id to read"
    arg_type = Int64
    default = 128
end
parsed_args = parse_args(tabargs)

const FILENAMES = parsed_args["filenames"]
const FILEIDMIN = parsed_args["fileidmin"]
const FILEIDMAX = parsed_args["fileidmax"]
@assert(FILEIDMIN <= FILEIDMAX)

const NFILES    = FILEIDMAX - FILEIDMIN + 1


##################################################
# Creating the sumup file
##################################################
const SUMUPFILENAME = FILENAMES*string(FILEIDMIN)*"_to_"*string(FILEIDMAX)*"_sumup.hf5"
sumupfile = h5open(SUMUPFILENAME,"w")
close(sumupfile)

##################################################
# Including IO functions for core 
##################################################
include("./IO.jl")

########################
# Return the file id given its position in 
# the list j in 1...NFILES
########################
function fileid(j::Int64)
    return FILEIDMIN - 1 + j 
end

########################
# Sumup function
#
# To specify before use :
#
# - otherattributes() : Return a list of specific attributes shared by all files to read in the first one
#
# - othersumup!(att_1stfile,listNTRY,TOTNTRY) :
#       Treat all the other specific things to sum up in the file (+ write in the sumup file)
#       
########################
function Sumup!()
    ##################################################
    # Quantities that are the same for every try/proc
    # + Only "seed" parameter is different in 
    # the files --> to sum up in an array and 
    # replace
    # + "ntry" should be replaced by the total number 
    # of try also
    ##################################################
    #####
    # General attributes (every simulation)
    #####
    att_toread = copy(att_list)
    deleteat!(att_toread, findall(x->x=="seed",att_toread)) # random seed is not constant over files
    deleteat!(att_toread, findall(x->x=="ntry",att_toread)) # we will save as "ntry" the total number of tries
    #####
    # Add other general attributes of interest (same for all files)
    #####
    otheratt = otherattributes()
    att_toread = vcat(att_toread,otheratt) # Add other general attributes to the list
    #####
    # Copy the general attributes from the first file
    #####
    att_1stfile = att_read(FILENAMES*string(fileid(1))*".hf5",att_toread)
    att_dump_fromdict!(SUMUPFILENAME,"r+",att_1stfile)

    #####
    # Saved time lists (to copy from 1 try in 1 file)
    #####
    firstfile = h5open(FILENAMES*string(fileid(1))*".hf5","r")
    ensaved_times = read(firstfile,"test_1_ensaved_t")
    testsaved_times = read(firstfile,"test_1_testsaved_t")
    fullsaved_times = read(firstfile,"test_1_fullsaved_t")
    close(firstfile)

    #####
    # Usefull attributes (arrays lenghts)
    #####
    NENSAVES   = att_1stfile["nensaves"]

    ############
    # Retrieving the total number of tries for memory allocation of variables of interest
    # And summarizing the random seed used in one list 
    ############
    NTRYlist = zeros(Int64,NFILES)
    SEEDlist = zeros(Int64,NFILES)
    for j = 1:NFILES
        file = h5open(FILENAMES*string(fileid(j))*".hf5","r")
        NTRYlist[j] = read(file,"ntry")
        SEEDlist[j] = read(file,"seed")
        close(file)
    end
    TOTNTRY = sum(NTRYlist)

    ##########
    # Total energy/momentum conservation 
    ##########
    Etotvar = zeros(Float64,NENSAVES,TOTNTRY)
    Etotinit= zeros(Float64,TOTNTRY) # Initial energy of each try
    Ptotvar = zeros(Float64,NENSAVES,TOTNTRY)
    Ptotinit= zeros(Float64,TOTNTRY) # Initial momentum of each try (in the mean DF reference frame not the simulation recentered one)
    ##########
    # Momentum center 
    ##########
    vmflist = zeros(Float64,TOTNTRY)

    trycount = 1
    for j = 1:NFILES
        file = h5open(FILENAMES*string(fileid(j))*".hf5","r")
        #####
        # Adding the i-th try of the j-th proc tabE and tabP to the sumup tab
        #####
        for i = 1:NTRYlist[j]
            #####
            # Conservation of the total energy/momentum
            #####
            actEtot = read(file,"test_"*string(i)*"_tabE")
            actPtot = read(file,"test_"*string(i)*"_tabP")
            Ein = actEtot[1]
            Etotinit[trycount] = Ein
            Pin = actPtot[1] 
            Ptotinit[trycount] = Pin
            for p = 1:NENSAVES
                Etotvar[p,trycount] = abs(actEtot[p]-Ein)
                Ptotvar[p,trycount] = abs(actPtot[p]-Pin)
            end
            ##########
            # Momentum center 
            ##########
            vmf = read(file,"test_"*string(i)*"_vmf")
            vmflist[trycount] = vmf
            #####
            trycount += 1
        end
        close(file)
    end
    #########
    # Dumping this concentrated infos in 1 file
    #########

    sumupfile = h5open(SUMUPFILENAME,"r+")

    write(sumupfile,"nfiles",NFILES)
    #########
    # Summed up informations
    #########
    write(sumupfile,"totntry",TOTNTRY)
    write(sumupfile,"ntrylist",NTRYlist)
    write(sumupfile,"seedlist",SEEDlist)

    write(sumupfile,"Etotvar",Etotvar)
    write(sumupfile,"Ptotvar",Ptotvar)
    write(sumupfile,"Etotinit",Etotinit)
    write(sumupfile,"Ptotinit",Ptotinit)

    write(sumupfile,"vmf",vmflist)

    write(sumupfile,"ensaved_t",ensaved_times)
    write(sumupfile,"testsaved_t",testsaved_times)
    write(sumupfile,"fullsaved_t",fullsaved_times)

    close(sumupfile)

    ##########
    # Specific sumup and dump in the sumup file
    ##########
    othersumup!(att_1stfile,NTRYlist,TOTNTRY)
end