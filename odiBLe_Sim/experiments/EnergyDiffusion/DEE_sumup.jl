include("./../../core/Sumup.jl")

##################################################
# Including DEE specific IO functions
##################################################
include("./DEE_IO.jl")

function otherattributes()
    return bins_params_list
end

function othersumup!(att::Dict{String,Any},
                    listntry::Array{Int64},
                    totntry::Int64)
    #####
    # Usefull attributes (arrays lenghts)
    #####
    NBINS = att["nbins"]
    NTIMESPERTRY = att["ntestsaves"]
    ###########
    # Quantities of interest : mean (\Delta E)^2 of each try for each bin
    # and bincount for each try/bin and mean energy in each bin/try
    # All tries sum up in only 1 array
    ###########
    meanDE2each = zeros(Float64,NBINS,NTIMESPERTRY,totntry)
    bincounteach = zeros(Int64,NBINS,totntry)
    meanEinbineach = zeros(Float64,NBINS,totntry)

    trycount = 1
    for j = 1:NFILES
        file = h5open(FILENAMES*string(fileid(j))*".hf5","r")
        for i = 1:listntry[j]
            #####
            # Adding the i-th try of the j-th proc mean (\Delta E)^2 to the sumup tab
            #####
            actmeanDE2      = read(file,"test_"*string(i)*"_meanDE2")
            actmeanEinbin   = read(file,"test_"*string(i)*"_meanEinbin")
            actbincount     = read(file,"test_"*string(i)*"_bincount")
            for q = 1:NTIMESPERTRY
                for p = 1:NBINS
                    meanDE2each[p,q,trycount[1]] = actmeanDE2[p,q]
                end
            end
            for p = 1:NBINS
                bincounteach[p,trycount[1]] = actbincount[p]
                meanEinbineach[p,trycount[1]] = actmeanEinbin[p]
            end
            #####
            trycount +=1
        end
        close(file)
    end

    #########
    # Dumping this concentrated infos in 1 file
    #########
    sumupfile = h5open(SUMUPFILENAME,"r+")
    write(sumupfile,"meanDE2",meanDE2each)
    write(sumupfile,"bincount",bincounteach)
    write(sumupfile,"meanEinbin",meanEinbineach)
    close(sumupfile)

end

Sumup!()