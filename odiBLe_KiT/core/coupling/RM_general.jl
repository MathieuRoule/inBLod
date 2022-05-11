if !(@isdefined structBasis_create)  # If it as never been included
    include("./Basis.jl")  # Basis handling structure and functions needed for response matrix computation
end
#####
const RM_RA_MAX	= BL_RA_MAX * Lambda # Maximal apocenter of integration for response matrix
##################################################
# Heaviside function on the interval [-1,1]
# Here, H[x] has a REAL argument, and returns
# 0   for x < -1
# 1/2 for x = -1
# 1   for -1 < x < 1
# 1/2 for x = 1
# 0   for 1 < x
# ATTENTION, the equality tests on Float64 might not be very robust
##################################################
function Heaviside(x::Float64)
    if     (x <  -1.0) # Left of the interval
        return 0.0
    elseif (x == -1.0) # Left edge of the interval
        return 0.5
    elseif (-1.0 < x < 1.0) # Within the interval
        return 1.0
    elseif (x ==  1.0) # Right edge of the interval
        return 0.5
    elseif (1.0 < x) # Right of the interval
        return 0.0
    end
end
######################################
# Computation of the response matrix
# using the splitting or the Legendre method 
# (common part of the algorithm)
######################################

# Choosing the method
if BLMETHOD == "Splitting"
    include("./RM_Splitting.jl") # To compute the response matrix using the splitting method
elseif BLMETHOD == "Legendre"
    include("./RM_Legendre.jl") # To compute the response matrix using the Legendre method
end


#################################################
# Computing the response matrix
# The cosine and the sine part of it
# by summing contribution of each frequency n in the FT
#################################################
function M_compute!(omg::Complex{Float64},structMat) 

    Mcos, Msin = structMat.Mcos, structMat.Msin

    fill!(Mcos,0.0 + 0.0*im) # Initialization of the matrix
    fill!(Msin,0.0 + 0.0*im) # Initialization of the matrix

    for n = (-NMAXRESMAT):NMAXRESMAT
        if ( mod(n,2) == 0 ) # if n is even, the FT of sine basis elmt is equal to 0 --> only cosine elmt are needed
            M_n!(n,omg,structMat,Mcos)
        else
            M_n!(n,omg,structMat,Msin)
        end
    end
end
