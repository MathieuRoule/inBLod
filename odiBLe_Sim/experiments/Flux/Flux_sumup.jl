include("./../../core/Sumup.jl")

##################################################
# Including Flux specific IO functions
##################################################
include("./../Binning/IO_E.jl")

function otherattributes()
    return bins_E_params_list
end

function othersumup!(att::Dict{String,Any},
                    listntry::Array{Int64},
                    totntry::Int64)
    NBINS_E = att["nbins_E"]
    NTIMESPERTRY = att["ntestsaves"]
    ###########
    # Quantities of interest : 
    # Total histogram of energy
    # All tries sum up in only 1 array
    ###########
    tothist_E = zeros(Int64,NBINS_E,NTIMESPERTRY)

    for j = 1:NFILES
        file = h5open(FILENAMES*string(fileid(j))*".hf5","r")
        #####
        # Quantities of interest
        #####
        acthist_E = read(file,"hist_E")
        #####
        # Adding all tries together
        #####
        for t = 1:NTIMESPERTRY
            for e = 1:NBINS_E
                tothist_E[e,t] += acthist_E[e,t]
            end
        end
        close(file)
    end

    #########
    # Dumping this concentrated infos in 1 file
    #########
    sumupfile = h5open(SUMUPFILENAME,"r+")
    write(sumupfile,"hist_E",tothist_E)
    close(sumupfile)

end

Sumup!()
