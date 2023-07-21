##################################################
# Function to save test particles states along simulation
# into tabtestsavex and tabtestsavev arrays
# we will dump all of them together at the end of the try (as we do for energy)
##################################################
function save_teststate!(k::Int64)
    testsaved_times[k]=TIME[1]
    for i = 1:NTESTPART
        tabtestsavex[i,k] = tabx[i]
        tabtestsavev[i,k] = tabv[i]
    end
end

##################################################
# At the end of each try 
# We dump (in addition with what is dumped 
# with dump_endtry from core IO):
# - The test states saved along the simulation
##################################################
function dump_endtry_testpart!(namefile::String,dataset::String,option::String="r+")
    #####
    file = h5open(namefile,option) # Opening the file with read and write
    #####
    write(file,dataset*"_testsaved_x",tabtestsavex) # Dumping the table of positions of test particles
    write(file,dataset*"_testsaved_v",tabtestsavev) # Dumping the table of positions of test particles
    #####
    close(file)
end


