# Loading all the needed packages
# and install them if they do not already exist

using Pkg # For packages management

const allpkg    =  [pkg.name for pkg in values(Pkg.dependencies())]
const pkglist   = ["ArgParse", "HDF5", "LinearAlgebra", "FastGaussQuadrature"]
const unknownpkg    = Array{String,1}()
for pkg in pkglist
    if !(pkg in allpkg)
        append!(unknownpkg,[pkg])
    end
end
if !(isempty(unknownpkg))
    println("The following packages are needed to use odiBLe_KiT code :")
    for pkg in unknownpkg
        println(pkg)
    end
    for pkg in unknownpkg
        Pkg.add(pkg)
    end
end



using ArgParse # To parse command-line arguments
using HDF5 # Writing results into HDF5 format files
using LinearAlgebra # Matrix inversion/determinant/eigenvalues for Response Matrix
using FastGaussQuadrature # Gauss Lengendre Quadrature