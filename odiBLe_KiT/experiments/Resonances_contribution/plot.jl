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
const nmax = read(file,"nmax")
const npart=read(file,"npart")

const tabDJJ = npart * read(file,"tabDJJ")

close(file)

gr()
pl=heatmap(1:nmax,1:nmax,log10.(tabDJJ),
            aspect_ratio = 1,
            colormap = :inferno, 
            clim = (-25,0), 
            colorbar_title = L"$\log (D_{JJ})$",
            interpolation = false)
title!(pl,state*"    "*coupling)
xlims!(0.5,nmax+0.5)
ylims!(0.5,nmax+0.5)
xlabel!(pl, L"$k\prime$")
ylabel!(pl, L"$k$")
savefig(pl,FILENAME*"_plot.png")