##################################################
# Compute the contribution 
# from the (k,kp) resonance to the 
# friction and diffusion coefficients and to the flux 
# given the orbits (ra,rap) and the frequency k*Omega(ra) (= kp*Omega(rap)).
# ATTENTION: We assume here that the frequency profile is 
# monotonous. Therefore (k,kp=k) resonances do not contribute to the flux
# as resonance condition is only fulfilled at J = Jp and then friction 
# and diffusion cancel each other.
# 3 outputs:
# (friction, diffusion, flux) contributions
##################################################
function get_resonance_contribution(k::Int64,kp::Int64,
                                    ra::Float64,rap::Float64,
                                    kdotomega::Float64,
                                    structMultipole=nothing,structBasis=nothing,structBL=nothing)
    #####
    # Prefactors 
    #####
    pref = 2.0*(pi)^(2)*mu # Prefactor for a given resonance (k,kp).
    prefdiff = 2.0 * pref  # ATTENTION, to the additional factor 2.0 specific to D_2
 
    #####
    # Computing the contribution of the (k,kp) resonance
    #####
    val_psiSQ = psiSQ(k,kp,ra,rap,kdotomega+0.0*im,structMultipole,structBasis,structBL)
    val_DF = DF(ra) # Value of DF in ra
    val_DFp = DF(rap) # Value of DF in rap
    val_dDFdJ = dDFdJ(ra) # Value of dF/dJ in ra
    val_dDFdJp = dDFdJ(rap) # Value of dF/dJ in rap
    val_dOmegadrap = dOmegadra(rap) # Value of dOmega/drap
    val_Jac = JacJtora(rap) # Value of the Jacobian |dJ/drap|
    #####
    fric = pref * (k*kp*val_Jac*val_psiSQ*val_dDFdJp)/(abs(kp*val_dOmegadrap)) # Contribution from the resonance to friction
    diff = prefdiff * (k^(2)*val_Jac*val_psiSQ*val_DFp)/(abs(kp*val_dOmegadrap)) # Contribution from the resonance to diffusion 
    if (k == kp)
        flux = 0.0	# Assuming monotonuous frequency profile, (k,k) resonances does not contribute to the flux
    else
        flux = pref * (k*val_Jac*val_psiSQ)/(abs(kp*val_dOmegadrap))  *  (kp * val_dDFdJp * val_DF - k * val_dDFdJ * val_DFp)  # Contribution from the resonance to the flux
    end
    return (fric, diff, flux)
end


########################################################################
# Diffusion, Friction and Flux Coefficients
########################################################################

##################################################
# Compute the diffusion coefficient,
# D_2, the friction coefficient, D_1, and the total flux 
# for a given orbital location
# ATTENTION, since the resonance (k,kp) contributes
# the same than (-k,-kp), we only account for the
# resonances with k>0, and add a factor 2 in the end
# ATTENTION, we consider only resonances with k!=0 and kp!=0
# ATTENTION, since Omega(ra)>0, we can limit ourselves
# to searching only for kp>0
# ATTENTION, we use the fact that Psi[k,kp] = 0
# if k and kp do not have the same parity
##################################################
function get_AJ_DJJ_Flux(ra::Float64,
                structMultipole=nothing,
                structBasis=nothing,
                structBL=nothing) 
    
    resfric = 0.0 # Container for the friction coefficient
    resdiff = 0.0 # Container for the diffusion coefficient
    resflux = 0.0 # Container for the flux
    #####
    val_Omega = Omegafromra(ra) # Current value of the resonance frequency
    #####
    for k=1:NMAX # Loop over the k. ATTENTION, we limit ourselves to k>0
        kdotomega = k*val_Omega # Current resonance frequency
        #####
        # Loop over the kp
        # ATTENTION, we limit ourselves to kp>0
        # ATTENTION, we limit ourselves to (k-kp) even
        # @IMPROVE -- the starting point for kp could be written in a more readable form
        for kp=(mod(k+1,2)+1):2:NMAX # Loop over kp. ATTENTION, the `step' is 2
            rap = get_Resonance(kdotomega,kp) # Solving for the resonance condition
            #####
            if (rap != nothing) # We have found a solution to the resonance condition
                rescontrib = get_resonance_contribution(k,kp,ra,rap,kdotomega,structMultipole,structBasis,structBL) # Contribution from the current resonance
                resfric += rescontrib[1]
                resdiff += rescontrib[2]
                resflux += rescontrib[3]
            end
            #####
        end
    end
    #####
    resfric *= 2.0 # Multiplying by 2, as we only summed over k>0
    resdiff *= 2.0 # Multiplying by 2, as we only summed over k>0
    resflux *= 2.0 # Multiplying by 2, as we only summed over k>0
    #####
    return (resfric, resdiff, resflux) # Output
end

####################################################
# Compute the diffusion coefficient in energy E (instead of action J) for a given location.
# We use the fact that 
# Mean((Delta E)^2) = (dE/dJ)^2 * Mean((Delta J)^2)
#                   = Omega()^2 * Mean((Delta J)^2) 
# where Omega() is equivalently a function of the apocenter (ra), the energy (E) or the action (J).
# Then, we deduce that
# D_E,E = Omega^2 * D_J,J
####################################################
function get_DEE(ra::Float64,
                structMultipole=nothing,
                structBasis=nothing,
                structBL=nothing)
    return (Omegafromra(ra))^(2) * get_AJ_DJJ_Flux(ra,structMultipole,structBasis,structBL)[2]
end

###################################################
# Compute DEE given DJJ at apocenter ra
###################################################
function get_DEE_fromra_DJJ(ra::Float64,DJJ::Float64)
    return (Omegafromra(ra))^(2) * DJJ
end

###################################################
# Compute DEE given Omega and DJJ
###################################################
function get_DEE_fromOm_DJJ(Omega::Float64,DJJ::Float64)
    return (Omega)^(2) * DJJ
end

###################################################
# Compute the flux given the diffusion and friction coefficients.
# Using Flux(J) = A(J) * F0(J) - 0.5 * D(J) * dF0dJ(J).
###################################################
function get_FluxJ_from_ra_AJ_DJJ(ra::Float64,AJ::Float64,DJJ::Float64)
    val_DF = DF(ra) # Value of F in rap
    val_dDFdJ = dDFdJ(ra) # Value of dF/dJ in rap
    #####
    return AJ * val_DF - 0.5 * DJJ * val_dDFdJ # Output
end
