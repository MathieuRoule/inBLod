if !(@isdefined M_compute!)
    include("./RM_general.jl") # Response matrix needed for the dressed prediction
end

# For substitution 
if !(@isdefined psi_Multipole)  # If it as never been included
    include("./Psi_Multipole.jl")  # Basis handling structure and functions
end

##################################################
# Dressed coupling coefficient
# using the basis and the splitting or Legendre method
# ATTENTION, this function returns 0 if (n-np) odd
# ATTENTION - If k is even, only the cosine elmt have non-zero FT value
# ATTENTION - If k is odd, only the sine elmt have non-zero FT value
# ATTENTION, if ra or rap > L, the interaction potential U
# is ill-represented by the basis. Therefore we switch to the 
# Multipole coupling coefficients.
##################################################
const IMat = zeros(Complex{Float64},NBASISELMT,NBASISELMT) # Static container for the identity matrix
##########
for p=1:NBASISELMT # Loop over the radial elements to fill in the identity matrix
    IMat[p,p] = 1.0 + 0.0im # Creating the identity matrix
end

const BL_FIRSTCALL = ones(Int64,1)

function psid_BL(k::Int64,kp::Int64,
               ra::Float64,rap::Float64,
               omgres::Complex{Float64}=0.0+0.0*im, 
               structMultipole=nothing, # For substitution
               structBasis=nothing, 
               structBL=nothing) 

    #####
    # If k and kp do not have the same parity
    # the coupling coefficient is equal to 0
    if (mod(k-kp,2) != 0)
        return 0.0
    end
    #####
    # If ra or rap > L, we switch to th Multipole
    # coupling coefficient
    if ((ra > L_per) || (rap > L_per))
        return psi_Multipole(k,kp,ra,rap,omgres,structMultipole,nothing)
    end
    
    #####
    if BL_FIRSTCALL[1] == 1
        BL_fullinit!(structBasis,structBL)
        BL_FIRSTCALL[1] = 0
    end

    #####
    # Giving a name to the arrays locally used
    tabpsicos, tabpsisin = structBasis.tabpsicos, structBasis.tabpsisin # Tables of the basis elmt value
    tabftpsicos, tabftpsisin = structBasis.tabftpsicos, structBasis.tabftpsisin # Tables of the Fourier transform basis elmt value for n
    tabftpsicosp, tabftpsisinp = structBasis.tabftpsicosp, structBasis.tabftpsisinp # Tables of the Fourier transform basis elmt value for np

    fourierbasisvalue!(k,ra,tabpsicos,tabpsisin,tabftpsicos,tabftpsisin) # n,ra fourier transform calculation
    fourierbasisvalue!(kp,rap,tabpsicos,tabpsisin,tabftpsicosp,tabftpsisinp) # np,rap fourier transform calculation

    #####
    # Part associated to the response matrix
    #####
    #= @timeit to "M_compute" =# M_compute!(omgres,structBL)
    
    if mod(k,2) == 0 
        ResMat = structBL.Mcos
        tabftpsi, tabftpsip = tabftpsicos, tabftpsicosp
    else
        ResMat = structBL.Msin
        tabftpsi, tabftpsip = tabftpsisin, tabftpsisinp
    end

    #= @timeit to "inversing" =# invloc = inv(Symmetric(IMat - ResMat)) # @ TO IMPROVE -- in place

    psid = 0.0 + 0.0*im
    for q = 1:NBASISELMT
        for p = 1:NBASISELMT
            psid -= tabftpsi[p]*invloc[p,q]*tabftpsip[q]
        end
    end

    return psid # Returning the coupling coefficient 
end
