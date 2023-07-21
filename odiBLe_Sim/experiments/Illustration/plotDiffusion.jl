###########################
# Basic program to plot 
# the diffusion illustration
# computed using RunDiffusion.jl and saved 
# in a HDF5 file.
# Name of the file to give in argument 
# when calling the code.
###########################

# Testing if all needed packages are installed 
# If not install them.
include("./../../pkgcheck.jl")
pkgcheck("ArgParse","HDF5","Plots","LaTeXStrings","ColorSchemes")

# Charging the needed packages
using ArgParse
using HDF5
using Plots
using LaTeXStrings
using ColorSchemes
using Base.Order # For array sorting purpose


# Reading file name parsed in command line
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--filename"
    help = "Path + name of the HDF5 file (without the extension)"
    arg_type = String
    default = "./KinTh"
end
parsed_args = parse_args(tabargs)
const FILENAME = parsed_args["filename"]

file=h5open(FILENAME*".hf5","r")
const state = read(file,"fbathsampling")
const coupling = read(file,"coupling")

const npart=read(file,"nbathpart")

const tabt = read(file,"test_1_testsaved_t")
const tabx = read(file,"test_1_testsaved_x")
const tabv = read(file,"test_1_testsaved_v")

close(file)

const ntestpart, ntimes = size(tabx)
const tabE = zeros(Float64,ntestpart,ntimes)

#####
# Define the potential
#####
include("./../../core/Constants.jl")
if state == "GTE"
    include("./../../core/states/GTE.jl")
    const pot = potential_GTE
elseif state == "Plummer"
    include("./../../core/states/Plummer.jl")
    const pot = potential_pol
elseif state == "Polytrope"
    include("./../../core/states/Polytrope.jl")
    const pot = potential_plum
else
    error("Unknown state")
end

for it = 1:ntimes
    for kp = 1:ntestpart
        v = tabv[kp,it]
        px = pot(tabx[kp,it])
        tabE[kp,it] = 0.5*(v^2) + px
    end
end 

const finalenergy_order = [k for k = 1:ntestpart]

sort!(finalenergy_order,QuickSort,Perm(ForwardOrdering(),tabE[:,end]))


plE=plot(tabt,transpose(tabE)[:,finalenergy_order],legend=false,lw = 2,palette=palette(:viridis,ntestpart))
xlabel!(plE, L"$t / t_{\mathrm{dyn}}$")
ylabel!(plE, L"$E / E_0$")

plxv=plot(transpose(tabx)[:,finalenergy_order],transpose(tabv)[:,finalenergy_order],legend=false,lw=0.5,palette=palette(:viridis,ntestpart))
xlabel!(plxv, L"$x / \Lambda$")
ylabel!(plxv, L"$v / \sigma$")

gatherplot = plot(plE, plxv,layout=(1,2),size=(900,400),plot_title="Diffusion   "*state*"   "*coupling,margin=4Plots.mm)
savefig(gatherplot,FILENAME*"_diffusion.pdf")