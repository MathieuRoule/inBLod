#####
# Work only with one sumup file
#####

using ArgParse # To be able to parse command-line arguments
using HDF5 # To be able to dump/read HDF5 files

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--filename"
    help = "Sumup file to read (path + name)"
    arg_type = String
end
parsed_args = parse_args(tabargs)

const FILENAME = parsed_args["filename"]

##################################################
# Including IO functions for core and DEE specific
##################################################
include("./../../core/IO.jl")
include("./DEE_IO.jl")

##################################################
# Creating the merged sumup file
##################################################
const AVGFILENAME = FILENAME*"_averaged.hf5"
avgfile = h5open(AVGFILENAME,"w")
close(avgfile)

#################################################
# Reading and copying the attributes
# dataset that don't need to be averaged (simple copy)
#################################################
const att_toread = copy(att_list)
deleteat!(att_toread, findall(x->x=="seed",att_toread)) # random seed is not constant over files
deleteat!(att_toread, findall(x->x=="ntry",att_toread)) # we will save as "ntry" the total number of tries
push!(att_toread,"ensaved_t")
push!(att_toread,"testsaved_t")
for binsatt in bins_params_list # Attributes linked to binning
    push!(att_toread,binsatt)
end
push!(att_toread,"nfiles")
push!(att_toread,"totntry")
push!(att_toread,"ntrylist")
push!(att_toread,"seedlist")

const file_att = att_read(FILENAME*".hf5",att_toread)

att_dump_fromdict!(AVGFILENAME,"r+",file_att)

#####
# Usefull attributes (arrays lenghts)
#####
const TOTNTRY = file_att["totntry"]
const NBINS = file_att["nbins"]
const NTIMESPERTRY = file_att["ntestsaves"]
const NENTOT = file_att["nensaves"]


#################################################
# Reading information from the sumed up file
# that need to be average
#################################################
file = h5open(FILENAME*".hf5","r")
#####
# DEE informations
#####
const meanDE2each = read(file,"meanDE2")
const bincounteach = read(file,"bincount")
const meanEinbineach = read(file,"meanEinbin")
#####
# Energy conservation
#####
const Etotvareach = read(file,"Etotvar")
const Ptotvareach = read(file,"Ptotvar")
const Etotinit = read(file,"Etotinit")
#####
close(file)





#################################################
# Averaging quantities over realizations
#################################################

#####
# Total mean of the weighted mean estimator calculation
# with the weight of 1 realization linked to the number of particles in the bin used to calculate the
# mean over the bin :
# w_k = bincount(bin,k) / ( sum_k (bincount(bin,k)) )
# Therefore, the mean is 
#
# m(bin) = sum_k (w_k * m_k)
# where m_k is the mean of DE^2 over the bin of the k-th realization.
#####
const totbincount = zeros(Int64,NBINS)
for r = 1:TOTNTRY
    for p = 1:NBINS
        totbincount[p] += bincounteach[p,r]
    end
end
const weights = zeros(Float64,NBINS,TOTNTRY)
for r = 1:TOTNTRY
    for p = 1:NBINS
        weights[p,r] = (totbincount[p] == 0) ? 0.0 : bincounteach[p,r] / totbincount[p]
    end
end

const meanDE2tot = zeros(Float64,NBINS,NTIMESPERTRY)
for r = 1:TOTNTRY
    for q = 1:NTIMESPERTRY
        for p = 1:NBINS
            meanDE2tot[p,q] += weights[p,r] * meanDE2each[p,q,r] 
        end
    end
end

const meanEinbintot = zeros(Float64,NBINS)
for r = 1:TOTNTRY
    for p = 1:NBINS
        meanEinbintot[p] += weights[p,r] * meanEinbineach[p,r]
    end
end

#####
# Unbiaised variance of a weighted sampling : (hyp : sum (w_k) = 1)
#
# s^2 = 1 / (1 - sum(w_k^2)) * sum (w_k * (m_k - m)^2)
# where m_k is the mean (on the bin not on realizations !) DE^2 of the k-th realization.
#####
const totweightsquared = zeros(Float64,NBINS)
for r = 1:TOTNTRY
    for p = 1:NBINS
        totweightsquared[p] += (weights[p,r])^(2)
    end
end
const s2pref = zeros(Float64,NBINS)
for p = 1:NBINS
    s2pref[p] = (totweightsquared[p] == 1.0) ? 0.0 : 1.0 / (1.0 - totweightsquared[p])
end

const varDE2tot = zeros(Float64,NBINS,NTIMESPERTRY)
for r = 1:TOTNTRY
    for q = 1:NTIMESPERTRY
        for p = 1:NBINS
            varDE2tot[p,q] += s2pref[p] * weights[p,r] * (meanDE2each[p,q,r] - meanDE2tot[p,q])^(2)
        end
    end
end

#####
# Error on the mean :
# SEM^2  =  s^2 / NTRY
# SEM = sqrt( varDE2 / NTRY )
#####
const SEMDE2tot = zeros(Float64,NBINS,NTIMESPERTRY)
for q = 1:NTIMESPERTRY
    for p = 1:NBINS
        SEMDE2tot[p,q] = sqrt(varDE2tot[p,q]/TOTNTRY)
    end
end 

#####
# Total energy/momentum conservation
####
const meanEtotvar = zeros(Float64,NENTOT)
const stdEtotvar = zeros(Float64,NENTOT)
const meanPtotvar = zeros(Float64,NENTOT)
const stdPtotvar = zeros(Float64,NENTOT)

for r = 1:TOTNTRY
    for p = 1:NENTOT
        meanEtotvar[p] += Etotvareach[p,r]
        meanPtotvar[p] += Ptotvareach[p,r]
    end
end
for p = 1:NENTOT
    meanEtotvar[p] /= TOTNTRY
    meanPtotvar[p] /= TOTNTRY
end
for r = 1:TOTNTRY
    for p = 1:NENTOT
        stdEtotvar[p] += (Etotvareach[p,r] - meanEtotvar[p])^2
        stdPtotvar[p] += (Ptotvareach[p,r] - meanPtotvar[p])^2
    end
end
for p = 1:NENTOT
    stdEtotvar[p] /= (TOTNTRY == 1) ? 1.0 : TOTNTRY - 1.0 # unbiaised variance
    stdPtotvar[p] /= (TOTNTRY == 1) ? 1.0 : TOTNTRY - 1.0 # unbiaised variance
end 
for p = 1:NENTOT
    stdEtotvar[p] = sqrt(stdEtotvar[p]) # standard deviation
    stdPtotvar[p] = sqrt(stdPtotvar[p]) # standard deviation
end 



#################################################
# Dumping the averaged information and standard 
# deviations in a file
#################################################

avgfile = h5open(AVGFILENAME,"r+")

#####
# Averaged informations
#####
write(avgfile,"meanDE2tot",meanDE2tot)
write(avgfile,"varDE2tot",varDE2tot)
write(avgfile,"SEMDE2tot",SEMDE2tot)
write(avgfile,"bincounttot",totbincount)
write(avgfile,"meanEinbintot",meanEinbintot)

##########
# Averaged Conservation of total energy/momentum
##########
write(avgfile,"meanEtotvar",meanEtotvar)
write(avgfile,"meanPtotvar",meanPtotvar)
write(avgfile,"stdEtotvar",stdEtotvar)
write(avgfile,"stdPtotvar",stdPtotvar)
write(avgfile,"Etotinit",Etotinit)

close(avgfile)

