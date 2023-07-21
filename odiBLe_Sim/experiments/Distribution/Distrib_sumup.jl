include("./../../core/Sumup.jl")

##################################################
# Including Distrib specific IO functions
##################################################
include("./IO_part.jl")

function otherattributes()
    return Array{String,1}()
end

function othersumup!(att::Dict{String,Any},
                    listntry::Array{Int64},
                    totntry::Int64)
    NPART = att["npart"]
    NTESTPART = att["ntestpart"]
    NTESTSAVES = att["ntestsaves"]
    NFULLSAVES = att["nfullsaves"]
    ###########
    # Quantities of interest : 
    # - testsaved x and v
    # - fullsaved x and v
    ###########
    tabtestx = zeros(Float64,NTESTPART,NTESTSAVES,totntry)
    tabtestv = zeros(Float64,NTESTPART,NTESTSAVES,totntry)
    tabfullx = zeros(Float64,NPART,NFULLSAVES,totntry)
    tabfullv = zeros(Float64,NPART,NFULLSAVES,totntry)

    trycount = 1
    for j = 1:NFILES
        file = h5open(FILENAMES*string(fileid(j))*".hf5","r")
        for i = 1:listntry[j]
            #####
            # Adding the i-th try of the j-th proc to the sumup tab test particles saved
            #####
            acttestx = read(file,"test_"*string(i)*"_testsaved_x")
            acttestv = read(file,"test_"*string(i)*"_testsaved_v")
            for q = 1:NTESTSAVES
                for p = 1:NTESTPART
                    tabtestx[p,q,trycount] = acttestx[p,q]
                    tabtestv[p,q,trycount] = acttestv[p,q]
                end
            end
            #####
            # Adding the i-th try of the j-th proc to the sumup tab full particles saved
            #####
            for q = 1:NFULLSAVES
                actfullx = read(file,"test_"*string(i)*"_t_"*string(q)*"_fullsaved_x")
                actfullv = read(file,"test_"*string(i)*"_t_"*string(q)*"_fullsaved_v")
                for p = 1:NPART
                    tabfullx[p,q,trycount] = actfullx[p]
                    tabfullv[p,q,trycount] = actfullv[p]
                end
            end
            trycount +=1
        end
        close(file)
    end

    #########
    # Dumping this concentrated infos in 1 file
    #########
    sumupfile = h5open(SUMUPFILENAME,"r+")
    write(sumupfile,"testsaved_x",tabtestx)
    write(sumupfile,"testsaved_v",tabtestv)
    write(sumupfile,"fullsaved_x",tabfullx)
    write(sumupfile,"fullsaved_v",tabfullv)
    close(sumupfile)

end

Sumup!()