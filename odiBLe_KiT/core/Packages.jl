# Loading all the needed packages
# and install them if they are not

include("./../pkgcheck.jl")
pkgcheck("ArgParse","HDF5","LinearAlgebra","FastGaussQuadrature")

using ArgParse # To parse command-line arguments
using HDF5 # Writing results into HDF5 format files
using LinearAlgebra # Matrix inversion/determinant/eigenvalues for Response Matrix
using FastGaussQuadrature # Gauss Lengendre Quadrature