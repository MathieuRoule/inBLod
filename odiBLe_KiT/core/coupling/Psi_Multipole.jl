if !(@isdefined structMultipole_create)
    include("./Multipole.jl") # Multipole method needed for the prediction
end

##################################################
# Bare coupling coefficient
# ATTENTION, this function returns 0 if (k-kp) odd
##################################################
function psi_Multipole(k::Int64,kp::Int64,
               ra::Float64,rap::Float64,
               omgres::Complex{Float64}=0.0+0.0*im, # No utility here
               structMultipole=nothing,
               structBasis=nothing, # No utility here
               structBL=nothing) # No utility here
    #####
    # If k and kp do not have the same parity
    # the coupling coefficient is equal to 0
    if (mod(k-kp,2) != 0)
        return 0.0
    end
    #####
    # Giving a name to the arrays locally used
    tabx, tabxp = structMultipole.tabx, structMultipole.tabxp # Tables of the positions
    tabg, tabgp = structMultipole.tabg, structMultipole.tabgp # Tables of the prefactors
    tabdeltaP0, tabdeltaP1 = structMultipole.tabdeltaP0, structMultipole.tabdeltaP1 # Tables for deltaP0 and deltaP1
    tabdeltaQ0, tabdeltaQ1 = structMultipole.tabdeltaQ0, structMultipole.tabdeltaQ1 # Tables for deltaQ0 and deltaQ1
    tabw = structMultipole.tabw # Table of the array sorting
    #####
    tabxtabg!(k,ra,tabx,tabg) # Updating positions and prefactors for the first particle
    tabxtabg!(kp,rap,tabxp,tabgp) # Updating positions and prefactors for the second particle
    tabw!(tabw,tabx,tabxp) # Getting the respective order of the particles
    #####
    # Computing the partial sums
    for j=1:K_FT
        #####
        # Computing tabdeltaP0 and tabdeltaP1
        deltaP0, deltaP1 = 0.0, 0.0 # Initialisation
        for i=(tabw[j]+1):(tabw[j+1]) # ATTENTION, indices are shifted by one
            g = tabg[i]
            deltaP0 += g
            deltaP1 += g*tabx[i]
        end
        #####
        tabdeltaP0[j] = deltaP0
        tabdeltaP1[j] = deltaP1
        #####
        # Computing tabdeltaQ0 and tabdeltaQ1
        deltaQ0, deltaQ1 = 0.0, 0.0
        for i=(tabw[j+1]+1):(tabw[j+2]) # ATTENTION, indices are shifted by one
            g = tabg[i]
            deltaQ0 += g
            deltaQ1 += g*tabx[i]
        end
        #####
        tabdeltaQ0[j] = deltaQ0
        tabdeltaQ1[j] = deltaQ1
    end
    #####
    # Computing the contributions from the sum P
    x0, x1 = tabdeltaP0[1], tabdeltaP1[1] # Initial values of deltaP0 and deltaP1
    sumP = tabgp[1]*(tabxp[1]*x0 - x1) # Initial contribution to P from j=1
    #####
    for j=1:(K_FT-1) # Loop over the samples
        x0 += tabdeltaP0[j+1] # Computing P0[j+1]
        x1 += tabdeltaP1[j+1] # Computing P1[j+1]
        #####
        sumP += tabgp[j+1]*(tabxp[j+1]*x0 - x1) # Contribution from j+1 to P
    end
    #####
    # Computing the contributions from the sum Q
    x0, x1 = tabdeltaQ0[K_FT], tabdeltaQ1[K_FT] # Initial values of deltaQ0 and deltaQ1
    sumQ = tabgp[K_FT]*(-tabxp[K_FT]*x0 + x1) # Initial contribution to P from j=K_FT
    #####
    for j=K_FT:-1:2 # Loop over the samples. ATTENTION, this is a decreasing recurrence
        x0 += tabdeltaQ0[j-1] # Computing Q0[j+1]
        x1 += tabdeltaQ1[j-1] # Computing Q1[j+1]
        #####
        sumQ += tabgp[j-1]*(-tabxp[j-1]*x0 + x1) # Contribution from j-1 to Q
    end
    #####
    psi = ((sumP+sumQ)*4.0*G)/(pi^(2)*K_FT^(2)) # Summing the result from P and Q, and rescaling by the correct prefactor
    return psi # Returning the coupling coefficient 
end
