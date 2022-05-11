if !(@isdefined structBasis_create)  # If it as never been included
    include("./Basis.jl")  # Basis handling structure and functions
end

# For substitution 
if !(@isdefined psi_Multipole)  # If it as never been included
    include("./Psi_Multipole.jl")  # Basis handling structure and functions
end

##################################################
# Non-dressed Bare coupling coefficient
# using the basis 
# ATTENTION, this function returns 0 if (k-kp) odd
# ATTENTION, if ra or rap > L, the interaction potential U
# is ill-represented by the basis. Therefore we switch to the 
# Multipole coupling coefficients.
##################################################
function psi_Basis_Landau(k::Int64,kp::Int64,
               ra::Float64,rap::Float64,
               omgres::Complex{Float64}=0.0+0.0*im, # No utility here
               structMultipole=nothing, # For substitution if ra or rap > L_per
               structBasis=nothing, 
               structBL=nothing) # No utility here
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
        return psi_Multipole(k,kp,ra,rap,omgres,structMultipole,nothing,nothing)
    end
    #####
    # Giving a name to the arrays locally used
    tabpsicos, tabpsisin = structBasis.tabpsicos, structBasis.tabpsisin # Tables of the basis elmt value
    tabftpsicos, tabftpsisin = structBasis.tabftpsicos, structBasis.tabftpsisin # Tables of the Fourier transform basis elmt value for k
    tabftpsicosp, tabftpsisinp = structBasis.tabftpsicosp, structBasis.tabftpsisinp # Tables of the Fourier transform basis elmt value for kp

    fourierbasisvalue!(k,ra,tabpsicos,tabpsisin,tabftpsicos,tabftpsisin) # n,ra fourier transform calculation
    fourierbasisvalue!(kp,rap,tabpsicos,tabpsisin,tabftpsicosp,tabftpsisinp) # np,rap fourier transform calculation

    if mod(k,2) == 0 
        tabftpsi, tabftpsip = tabftpsicos, tabftpsicosp
    else
        tabftpsi, tabftpsip = tabftpsisin, tabftpsisinp
    end

    psi = 0.0
    for p = 1:NBASISELMT
        psi -= tabftpsi[p]*tabftpsip[p]
    end

    return psi # Returning the coupling coefficient
end
