##################################################
# Dressed coupling coefficient
# using the Balescu-Lenard, Landau-Multipole and Landau-Basis predictions 
# to complete the dressed part on the basis with the undressed part 
# on unrepresented smaller separations (using the difference between bare multipole and basis)
# PdS : "Pain de Sucre" 
# -> Psi^d_true = sum_{p=0}^{+\infty} sum_{q=0}^{+\infty} \psi^(p) (I - M)^{-1} \psi^(q)     (\psi^(p) basis elements)
#               = sum_{p=0}^{nbasiselt} sum_{q=0}^{nbasiselt} \psi^(p) (I - M)^{-1} \psi^(q) + sum_{unrepresented} \psi^(p) (I - M)^{-1} \psi^(q)
#               ~ Psi^d_BL  +  sum_{unrepresented} \psi^(p) I \psi^(q)                       (Collectif effects negligible at small scales)
#               = Psi^d_BL  + Psi_Multipole - Psi_Basis_Landau
#
# ATTENTION, if ra or rap > L, the interaction potential U
# is ill-represented by the basis. Therefore we switch to the 
# Multipole coupling coefficients.
##################################################
if !(@isdefined psi_Multipole)
    include("./Psi_Multipole.jl")
end
if !(@isdefined psi_Basis_Landau)
    include("./Psi_Basis.jl")
end
if !(@isdefined psid_BL)
    include("./Psid_BL.jl")
end

function psid_BL_PdS(k::Int64,kp::Int64,
               ra::Float64,rap::Float64,
               omgres::Complex{Float64}=0.0+0.0*im,
               structMultipole=nothing,
               structBasis=nothing, 
               structBL=nothing)
    
    #####
    # If ra or rap > L, we switch to th Multipole
    # coupling coefficient
    # (equivalent to what follow as it is also replaced in psi_BL and psi_Basis_Landau 
    # but prevent from calling psi_Multipole 3 times)
    if ((ra > L_per) || (rap > L_per))
        return psi_Multipole(k,kp,ra,rap,omgres,structMultipole,nothing,nothing)
    end

    psidBL = psid_BL(k,kp,ra,rap,omgres,structMultipole,structBasis,structBL)
    psiMulti = psi_Multipole(k,kp,ra,rap,omgres,structMultipole,nothing,nothing)
    psiBasis = psi_Basis_Landau(k,kp,ra,rap,omgres,structMultipole,structBasis,nothing)
    return psidBL + psiMulti - psiBasis # "Pain de Sucre"
end