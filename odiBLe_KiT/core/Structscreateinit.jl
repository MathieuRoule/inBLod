######
# Creation of the structures needed for computation
######
if COUPLING == "Landau-Multipole"
    const structMultipole = [structMultipole_create() for k=1:NTHREADS]
    const structBasis = [nothing for k=1:NTHREADS]
    const structBL = [nothing for k=1:NTHREADS]
elseif COUPLING == "Landau-Basis"
    const structMultipole = [structMultipole_create() for k=1:NTHREADS]
    const structBasis = [structBasis_create() for k=1:NTHREADS]
    const structBL = [nothing for k=1:NTHREADS]
elseif COUPLING == "Balescu-Lenard"
    const structMultipole = [structMultipole_create() for k=1:NTHREADS]
    const structBasis = [structBasis_create() for k=1:NTHREADS]
    if BLMETHOD == "Splitting"
        const structBL = [structSplitting_create() for k=1:NTHREADS]
    elseif BLMETHOD == "Legendre"
        const structBL = [structLegendre_create() for k=1:NTHREADS]
    end
elseif COUPLING == "BL-PdS"
    const structMultipole = [structMultipole_create() for k=1:NTHREADS]
    const structBasis = [structBasis_create() for k=1:NTHREADS]
    if BLMETHOD == "Splitting"
        const structBL = [structSplitting_create() for k=1:NTHREADS]
    elseif BLMETHOD == "Legendre"
        const structBL = [structLegendre_create() for k=1:NTHREADS]
    end
end

#####
#  Initialization for Balescu-Lenard (or BL-PdS)
#####
if (COUPLING == "Balescu-Lenard") || (COUPLING == "BL-PdS")
	BL_fullinit!(structBasis[1],structBL[1])
    if BLMETHOD == "Splitting"
        Threads.@threads for k=2:NTHREADS
            for j=1:K_BL
                structBL[k].tabomega[j] = structBL[1].tabomega[j] 
                structBL[k].tabdomegadra[j] = structBL[1].tabdomegadra[j] 
                for i=1:(2*NMAXRESMAT+1)
                    for q=1:NBASISELMT
                        for p=1:NBASISELMT
                            structBL[k].matG0[p,q,i,j] = structBL[1].matG0[p,q,i,j]
                        end
                    end
                end
            end 
        end
    elseif BLMETHOD == "Legendre"
        Threads.@threads for i=2:NTHREADS
            for om=1:4
                structBL[i].tabOmegapart[om] = structBL[1].tabOmegapart[om]
            end
            for k=1:K_BL
                for n=1:(2*NMAXRESMAT+1)
                    for q=1:NBASISELMT
                        for p=1:NBASISELMT
                            structBL[i].tabapqnk[p,q,n,k] = structBL[1].tabapqnk[p,q,n,k]
                        end
                    end
                end
            end
        end
    end
    BL_FIRSTCALL[1] = 0
end