######################################
# Computation of the response matrix
# using the splitting method 
# with the variable y = ra
######################################
struct structSplitting_type
    matG0 ::Array{Float64,4} # All the Matrix of g(y_i) the numerator in the integral (for every splitting point and every n)
    tabomega ::Array{Float64,1} # All the Omega(ra) for every splitting point
    tabdomegadra ::Array{Float64,1} # All the dOmegadra(ra) for every splitting point
    Mcos ::Array{Complex{Float64},2} # Response Matrix for cosine elmt (sum over n)
    Msin ::Array{Complex{Float64},2} # Response Matrix for sine elmt (sum over n)
end
##################################################
# Creates a RespMat_type
##################################################
function structSplitting_create()
    return structSplitting_type(zeros(Float64,(NBASISELMT,NBASISELMT,2*NMAXRESMAT+1,K_BL)), # g0 matrix at a specific point
                            zeros(Float64,K_BL),zeros(Float64,K_BL), # Omega(ra) and dOmegadra(ra) 
                            zeros(Complex{Float64},(NBASISELMT,NBASISELMT)),zeros(Complex{Float64},(NBASISELMT,NBASISELMT))) # cosine and sine matrix
end


##################################################
# Aleph_0 functions (different Im[omega] cases)
# Aleph_0(w) = INT[1/(u-w),{u,-1,1}],
# where the integral follows Landau's prescription
# As such, one has
# Aleph_0(w) = ln(1-w) - ln(-1-w)                   for Im[w] > 0
#            = ln|1-w| - ln|-1-w|  +   i*pi*Heaviside[w]      for Im[w] = 0
#            = ln(1-w) - ln(-1-w)  + 2*i*pi*Heaviside[Re[w]]  for Im[w] < 0
# Here, the three different cases correspond
# respectively to unstable; neutral; damped modes
##################################################
function aleph0_UNSTABLE(omg::Complex{Float64})
    return log(1.0-omg) - log(-1.0-omg)
end

function aleph0_NEUTRAL(omg::Complex{Float64})
    romg = real(omg)    
    return log(abs(1.0-romg)) - log(abs(-1.0-romg)) + Heaviside(romg)*pi*im
end

function aleph0_DAMPED(omg::Complex{Float64})
    romg = real(omg)    
    return log(1.0-omg) - log(-1.0-omg) + 2.0*Heaviside(romg)*pi*im
end
##################################################
# Choosing the function aleph0 depending
# on the imaginary part of the frequency
##################################################
function aleph0(omg::ComplexF64)
    # check the imaginary sign. if negative, use damped integration
    if (imag(omg) < 0.0)
        aleph0_DAMPED(omg)
    # if exactly zero, use neutral mode calculation
    elseif (imag(omg) == 0.0)
        aleph0_NEUTRAL(omg)
    # by default, use unstable integration
    else
        aleph0_UNSTABLE(omg)
    end
end

###################################################
# Difficult part of the response matrix computation.
# Computing the numerator value for every point ra
# and for every frequency in the FT n. All this part does not 
# depend on te resonance frequency omega.
# ATTENTION - If n is even, only the cosine elmt have non-zero FT value
# ATTENTION - If n is odd, only the sine elmt have non-zero FT value
###################################################
function Splitting_fullinit!(structBasis::structBasis_type,structBL)

    tabpcos, tabpsin = structBasis.tabpsicos, structBasis.tabpsisin # Tables of the basis elmt value
    tabftpcos, tabftpsin = structBasis.tabftpsicos, structBasis.tabftpsisin # Tables of the FT basis elmt value
    MatG0 = structBL.matG0
    tabomega, tabdomegadra = structBL.tabomega, structBL.tabdomegadra

    Delta_ra = RM_RA_MAX / K_BL # Integration step

    for istep = 1:K_BL 

        ra = (2.0 * istep - 1.0) * Delta_ra / 2.0 # Nodes of integration
        pref = - 2.0 * pi * JacJtora(ra) * dDFdJ(ra) # to be multiplied by n
        
        for n = -NMAXRESMAT:NMAXRESMAT # in matG0 the corresponding index is NMAXRESMAT+1 + n (since index are from 1 to 2*NMAXRESMAT+1 )
            if n == 0
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        MatG0[p,q,n+(NMAXRESMAT+1),istep] = 0.0
                    end
                end
            else
                #= @timeit to "fourier" =# fourierbasisvalue!(n,ra,tabpcos,tabpsin,tabftpcos,tabftpsin) # Computing the value of the FT of the basis elmt in x
                if mod(n,2) == 0 # if n is even, the FT of sine basis elmt is equal to 0 --> only cosine elmt are needed
                    for q = 1:NBASISELMT
                        for p = 1:NBASISELMT
                            MatG0[p,q,n+(NMAXRESMAT+1),istep] = n * pref * tabftpcos[p] * tabftpcos[q]
                        end
                    end
                else
                    for q = 1:NBASISELMT
                        for p = 1:NBASISELMT
                            MatG0[p,q,n+(NMAXRESMAT+1),istep] = n * pref * tabftpsin[p] * tabftpsin[q]
                        end
                    end
                end
            end
        end

        tabomega[istep] = Omegafromra(ra)
        tabdomegadra[istep] = dOmegadra(ra)
    end
end

##################################################
# Computing the part of response matrix associated
# to the n-th frequency in the FT
##################################################

function M_n_splitting!(n::Int64,omg::Complex{Float64},
                        structBL,
                        Mmat::Array{Complex{Float64},2}) 
    
    if (n == 0)
        return nothing
    end

    g0mat = structBL.matG0 
    tabomega, tabdomegadra = structBL.tabomega, structBL.tabdomegadra

    Delta_ra = RM_RA_MAX / K_BL # Integration step

    for istep = 1:K_BL 
        ra = (2.0 * istep - 1.0) * Delta_ra / 2.0 # Nodes of integration
        
        h0 = n*tabomega[istep]  # h0 = n*Omega(ra)
        h1 = n*tabdomegadra[istep]  # h1 = n*dOmegadra(ra)

        varpi = 2.0 * (omg - h0) / (Delta_ra * abs(h1)) # \varpi definition

        pref = aleph0(varpi) / abs(h1)

        #= @timeit to "mat" =# begin
            for q = 1:NBASISELMT
                for p = 1:NBASISELMT
                    Mmat[p,q] += pref * g0mat[p,q,n+(NMAXRESMAT+1),istep]
                end 
            end
        end
    end
end

##################################################
# Defining the response matrix method and initialization to Splitting
# ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
##################################################
const M_n! = M_n_splitting!
const BL_fullinit! = Splitting_fullinit!
