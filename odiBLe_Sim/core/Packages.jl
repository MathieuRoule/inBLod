# Loading all the needed packages
# and install them if they are not

# Testing if all needed packages are installed 
# If not install them.
include("./../pkgcheck.jl")
pkgcheck("ArgParse","HDF5","TimerOutputs")

using Base.Order # For array sorting purpose
using ArgParse # To be able to parse command-line arguments
using Random # To be able to fix the random seeds
using HDF5 # To be able to dump/read HDF5 files
using TimerOutputs # For algorithm benchmark
