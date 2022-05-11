##################################################
# Choice of susceptibility coefficient 
# way of calculation depending on coupling
##################################################

if COUPLING == "Landau-Multipole"
    include("./coupling/Psi_Multipole.jl") # Bare coupling coefficient using the multipole method
    const coupling_psi = psi_Multipole
elseif COUPLING == "Landau-Basis"
    include("./coupling/Psi_Basis.jl") # Bare coupling coefficient using the basis
    const coupling_psi = psi_Basis_Landau
elseif COUPLING == "Balescu-Lenard"
    include("./coupling/Psid_BL.jl") # Dressed coupling coefficient
    const coupling_psi = psid_BL
elseif COUPLING == "BL-PdS"
    include("./coupling/Psid_BL_PdS.jl") # Dressed coupling coefficient + using "Pain de Sucre" to smooth basis truncation (using Multipole for undressed contribution of truncated basis elements in the coupling coefficients)
    const coupling_psi = psid_BL_PdS
end

function psiSQ(n::Int64,np::Int64,
               ra::Float64,rap::Float64,
               omgres::Complex{Float64}=0.0+0.0*im,
               structMultipole=nothing,
               structBasis=nothing, 
               structBL=nothing)
    return (abs(coupling_psi(n,np,ra,rap,omgres,structMultipole,structBasis,structBL)))^(2)
end
