###########################
# Basic program to plot 
# the diffusion, friction and flux 
# coefficients computed using 
# RunBL.jl and saved 
# in a HDF5 file.
# Name of the file to give in argument 
# when calling the code.
###########################

# Testing if all needed packages are installed 
# If not install them.
include("./../../pkgcheck.jl")
pkgcheck("ArgParse","HDF5","Plots","LaTeXStrings")

# Charging the needed packages
using ArgParse
using HDF5
using Plots
using LaTeXStrings


# Reading file name parsed in command line
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--filename"
    help = "Path + name of the HDF5 file where the info to plot are"
    arg_type = String
    default = "./KinTh"
end
parsed_args = parse_args(tabargs)
const FILENAME = parsed_args["filename"]

file=h5open(FILENAME*".hf5","r")
const state = read(file,"state")
const coupling = read(file,"coupling")

const npart=read(file,"npart")

const tabJ = read(file,"tabJ")
const tabAJ = npart * read(file,"tabAJ")
const tabDJJ = npart * read(file,"tabDJJ")
const tabFlux = npart * read(file,"tabFlux")

const tabE = read(file,"tabE")
const tabDEE = npart * read(file,"tabDEE")

close(file)

pl=plot(tabJ, [tabAJ, tabDJJ, tabFlux], title = state*"    "*coupling, 
    label = ["Friction" "Diffusion" "Flux"], lw = 2)
xlabel!(pl, L"$J$")
ylabel!(pl, L"$N \times (A_J,\; D_{J,J},\; \mathcal{F}_J)$")
savefig(pl,FILENAME*"_plot.pdf")

plDEE=plot(tabE, tabDEE, title = state*"    "*coupling, 
    label = "Energy Diffusion", lw = 2, reuse = false)
xlabel!(plDEE, L"$E / E_0$")
ylabel!(plDEE, L"$N \times D_{EE}$")
savefig(plDEE,FILENAME*"_plotDEE.pdf")


