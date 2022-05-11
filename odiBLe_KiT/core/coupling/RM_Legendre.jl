######################################
# Computation of the response matrix
# using the Legendre method
######################################
struct structLegendre_type
    ######
    # Related to the considered points
    ######
    tabOmegapart ::Array{Float64,1} # Omega_0, Omega_L, Omega_m and Omega_e
    tabuGL ::Array{Float64,1}    # Table of the nodes   for the G-L quadrature
    tabwGL ::Array{Float64,1}    # Table of the weights for the G-L quadrature
    tabOmGLposn ::Array{Float64,1}    # Table of Omega on the nodes u_i for the G-L quadrature for resonances n > 0
    tabOmGLnegn ::Array{Float64,1}    # Table of Omega on the nodes u_i for the G-L quadrature for resonances n < 0
    tabraGLposn ::Array{Float64,1}    # Table of ra on the nodes u_i    for the G-L quadrature for resonances n > 0
    tabraGLnegn ::Array{Float64,1}    # Table of ra on the nodes u_i    for the G-L quadrature for resonances n < 0
    ######
    # Related to the Legendre polynomials
    ######
    tabINVcGL ::Array{Float64,1} # Table of the 1/c_k = (2k+1)/2 associated with the normalisation of the Legendre polynomials for the G-L quadrature
    tabPGL ::Array{Float64,2}    # Table of P[k,i]=P_k(u_i) used in the G-L quadrature.
    ######
    # Related to the coefficents of the response matrix
    ######
    tabGpqnk ::Array{Float64,4}  # Table of the value of G^pq_n(u_k) used for coefficient calculation
    tabapqnk ::Array{Float64,4}  # Table of the k-th coefficients of Legendre decomposition for every matrix element (p,q) and every resonnance (n)
    ######
    # Evaluation on the complex plane (integration with Landau's prescription of the Legendre polynomials)
    ######
    tabPLeg::Array{Complex{Float64},1} # Static container for tabPLeg
    tabQLeg::Array{Complex{Float64},1} # Static container for tabQLeg
    tabDLeg::Array{Complex{Float64},1} # Static container for tabDLeg
    ######
    # Response Matrix for cosine and sine basis elmts
    ######
    Mcos ::Array{Complex{Float64},2} # Response Matrix for cosine elmt (sum over n)
    Msin ::Array{Complex{Float64},2} # Response Matrix for sine elmt (sum over n)
end

##################################################
# Creates a Legendre_type
##################################################
function structLegendre_create()
    return structLegendre_type(zeros(Float64,4), # Omega0 for frequency normalization
                                zeros(Float64,K_BL),zeros(Float64,K_BL), # Tables of nodes u and weights w
                                zeros(Float64,K_BL),zeros(Float64,K_BL), # Tables of Omega on the nodes for n > 0 and n < 0
                                zeros(Float64,K_BL),zeros(Float64,K_BL), # Tables of ra on the nodes for n > 0 and n < 0
                                zeros(Float64,K_BL), # Tables of the 1/c_k
                                zeros(Float64,K_BL,K_BL), # Tables of the 1/c_k
                                zeros(Float64,NBASISELMT,NBASISELMT,2*NMAXRESMAT+1,K_BL), # Table of G (numerator) values
                                zeros(Float64,NBASISELMT,NBASISELMT,2*NMAXRESMAT+1,K_BL), # Table of the Legendre decomposition coefficients 
                                zeros(Complex{Float64},K_BL),zeros(Complex{Float64},K_BL),zeros(Complex{Float64},K_BL), # tabPLeg, tabQLeg and tabDLeg
                                zeros(Complex{Float64},(NBASISELMT,NBASISELMT)),zeros(Complex{Float64},(NBASISELMT,NBASISELMT))) # Mcos and Msin (Response Matrix)
end



##################################################
# Initialization of particular Omega :
# Omega_0 = Omega (0.0),
# Omega_L = Omega (J_L) the trucation frequency,
# Omega_m = (Omega_0 + Omega_L)/2 , the mean frequency and
# Omega_e =  (Omega_0 - Omega_L)/2 the extension frequency
##################################################
function Omegapart_init!(structLeg::structLegendre_type)
    Om0 = Omegafromra(0.0)
    OmL = Omegafromra(RM_RA_MAX)
    OmM = (Om0 + OmL) / 2.0
    OmE = (Om0 - OmL) / 2.0

    tabompart = structLeg.tabOmegapart
    tabompart[1] = Om0
    tabompart[2] = OmL
    tabompart[3] = OmM
    tabompart[4] = OmE
end
##################################################
# Function to initialize the nodes (u) and weights (w)
# of the Gauss-Legendre quadrature
# and then using the nodes : initializing Omega and ra on nodes (usefull for function computation).
# Indeed, we have
# u = sign(n) * (Omega - Omega_m) / Omega_e ,
# then
# Omega = Omega_m + sign(n) * Omega_e * u.
# We find the corresponding ra, using the bisection algorithm
# coded for resonances
##################################################
function GLpoints_init!(structLeg::structLegendre_type)

    # Nodes and weights
    tabuGL, tabwGL = structLeg.tabuGL, structLeg.tabwGL
    tabuGL[:], tabwGL[:] = gausslegendre(K_BL) # Computing the nodes (u) and weights (w) of the G-L quadrature

    # Corresponding Omega and ra
    tabompart = structLeg.tabOmegapart
    OmM, OmE = tabompart[3], tabompart[4]

    tabOmGLposn, tabOmGLnegn = structLeg.tabOmGLposn, structLeg.tabOmGLnegn
    tabraGLposn, tabraGLnegn = structLeg.tabraGLposn, structLeg.tabraGLnegn
    
    for k = 1:K_BL
        tabOmGLposn[k] = OmM + OmE * tabuGL[k]  # Omega(u)
        tabOmGLnegn[k] = OmM - OmE * tabuGL[k] # Omega(u)

        tabraGLposn[k] = get_Resonance(tabOmGLposn[k],1) # ra such that Omega(ra) = Omega(u)
        tabraGLnegn[k] = get_Resonance(tabOmGLnegn[k],1) # ra such that Omega(ra) = Omega(u)
    end
end
##################################################
# Function to initialise the normalisation constant
# INVc = 1/c_k = (2k+1)/2
# ATTENTION, it corresponds to the INVERSE.
##################################################
function GLinvc_init!(structLeg::structLegendre_type)
    tabINVcGL = structLeg.tabINVcGL
    for k=0:(K_BL-1) # Loop over the Legendre index
        tabINVcGL[k+1] = (2.0*k+1.0)/(2.0) # Filling in the value. ATTENTION, to the shift in the array index
    end
end
##################################################
# Function to initialise the values of the Legendre polynomials,
# as stored in tabPGLquad = P[k,i] = P_k(u_i)
##################################################
function GLPvalues_init!(structLeg::structLegendre_type)
    tabuGL, tabwGL = structLeg.tabuGL, structLeg.tabwGL
    tabPGL = structLeg.tabPGL
    for i=1:K_BL # Loop over the nodes
        u = tabuGL[i] # Current value of the node
        #####
        v0 = 1.0 # Initialisation of P_0(u)
        v1 = u   # Initialisation of P_1(u)
        #####
        tabPGL[1,i] = v0 # Filling in the value of P_0(u). ATTENTION, to the shift in the array index.
        tabPGL[2,i] = v1 # Filling in the value of P_1(u). ATTENTION, to the shift in the array index.
        #####
        for k=2:(K_BL-1) # Loop over the index 2 <= k
            v = (2.0*k-1.0)/(k)*u*v1 - (k-1.0)/(k)*v0 # Bonnet's recc. relation to get P_k(u)
            #####
            tabPGL[k+1,i] = v # Filling in the value P_k(u). ATTENTION, to the shift in the array index for k
            #####
            v0, v1 = v1, v # Shifting the temporary variables
        end
    end
end



###################################################
# Difficult part of the response matrix computation.
# Computing the values in each interpolation point (K_BL numbers)
# of the numerator for every matrix element (NBASISELMT*NBASISELMT)
# and for every resonance frequency (2*NMAXRESMAT+1).
# ATTENTION - If n is even, only the cosine elmt have non-zero FT value
# ATTENTION - If n is odd, only the sine elmt have non-zero FT value
###################################################
function Gpqnk_compute_init!(structBasis::structBasis_type,structLeg::structLegendre_type)

    # From structBasis
    tabpcos, tabpsin = structBasis.tabpsicos, structBasis.tabpsisin # Tables of the basis elmt value
    tabftpcos, tabftpsin = structBasis.tabftpsicos, structBasis.tabftpsisin # Tables of the FT basis elmt value
    
    # From structLeg
    tabu = structLeg.tabuGL
    tabraposn, tabranegn = structLeg.tabraGLposn, structLeg.tabraGLnegn
    Gpqnk = structLeg.tabGpqnk

    for k = 1:K_BL 

        u = tabu[k]

        # Negative resonances n < 0

        ra = tabranegn[k]
        pref = - 2.0 * pi * JacJtora(ra) * dDFdJ(ra) / dOmegadra(ra) # n < 0 --> sign(n) = -1       ## ATTENTION --> possible numerical divergence in ra -> 0 
        
        for n = -NMAXRESMAT:(-1) 
            #= @timeit to "fourier" =# fourierbasisvalue!(n,ra,tabpcos,tabpsin,tabftpcos,tabftpsin) # Computing the value of the FT of the basis elmt in x
            if mod(n,2) == 0 # if n is even, the FT of sine basis elmt is equal to 0 --> only cosine elmt are needed
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        Gpqnk[p,q,n+(NMAXRESMAT+1),k] = pref * tabftpcos[p] * tabftpcos[q]
                    end
                end
            else
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        Gpqnk[p,q,n+(NMAXRESMAT+1),k] = pref * tabftpsin[p] * tabftpsin[q]
                    end
                end
            end
        end

        # Positive resonances n > 0 

        ra = tabraposn[k]
        pref = 2.0 * pi * JacJtora(ra) * dDFdJ(ra) / dOmegadra(ra) # n > 0 --> sign(n) = +1

        for n = 1:NMAXRESMAT
            #= @timeit to "fourier" =# fourierbasisvalue!(n,ra,tabpcos,tabpsin,tabftpcos,tabftpsin) # Computing the value of the FT of the basis elmt in x
            if mod(n,2) == 0 # if n is even, the FT of sine basis elmt is equal to 0 --> only cosine elmt are needed
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        Gpqnk[p,q,n+(NMAXRESMAT+1),k] = pref * tabftpcos[p] * tabftpcos[q]
                    end
                end
            else
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        Gpqnk[p,q,n+(NMAXRESMAT+1),k] = pref * tabftpsin[p] * tabftpsin[q]
                    end
                end
            end
        end

        # n = 0
        for q = 1:NBASISELMT
            for p = 1:NBASISELMT
                Gpqnk[p,q,NMAXRESMAT+1,k] = 0.0
            end
        end
    end
end


###################################################.
# Computing the coordinates (K_BL numbers) of the Legendre polynomial decomposition of G^pq_n
# the numerator for every matrix element (NBASISELMT*NBASISELMT)
# and for every resonance frequency (2*NMAXRESMAT+1).
# ATTENTION - If n is even, only the cosine elmt have non-zero FT value
# ATTENTION - If n is odd, only the sine elmt have non-zero FT value
###################################################
function apqnk_compute_init!(structLeg::structLegendre_type)
    
    # From structLeg
    tabu, tabw = structLeg.tabuGL, structLeg.tabwGL # Points and weights
    tabinvc = structLeg.tabINVcGL
    PGL = structLeg.tabPGL # Values of the Legendre polynomial on the interpolation points  P[k,i] = P_k(u_i)
    apqnk, Gpqnk = structLeg.tabapqnk, structLeg.tabGpqnk # Coefficients of the Legendre decompostion and values of the numerator on every interpolation point

    fill!(apqnk,0.0) # Initialization

    for k = 1:K_BL # Over the coordinates (equivalently over the Legendre polynomials)

        invck = tabinvc[k]

        for i = 1:K_BL # Over the interpolation points

            wi = tabw[i]
            Pki = PGL[k,i]

            for n = 1:(2*NMAXRESMAT+1) # over the resonances
                for q = 1:NBASISELMT
                    for p = 1:NBASISELMT
                        apqnk[p,q,n,k] += invck * wi * Pki * Gpqnk[p,q,n,i]
                    end
                end
            end
        end

    end
end


###################################################
# Full initialization
###################################################
function Legendre_fullinit!(structBasis::structBasis_type,structBL)
    # Particular omega initialization
    Omegapart_init!(structBL)
    # Points initialization
    GLpoints_init!(structBL) 
    # Inverse of normalization coefficient initialization
    GLinvc_init!(structBL)
    # P[k,i] = P_k(u_i) values of the Legendre polynomials on the considered points
    GLPvalues_init!(structBL)
    # Computing the value of the numerator G^pq_n(u_k) on the interpolation points
    Gpqnk_compute_init!(structBasis,structBL)
    # Computing the value of the numerator Legendre components a^pq_n(k) (to be * by P_k(u))
    apqnk_compute_init!(structBL)
end


##################################################
# Function to fill in all the Legendre arrays
# for an UNSTABLE mode, i.e. Im[w] > 0.0
##################################################
function tabLeg!_UNSTABLE(omg::Complex{Float64},
                          struct_tabLeg::structLegendre_type)
    ##########
    tabDLeg = struct_tabLeg.tabDLeg # Name of the container for D_k(w)
    ##########
    # Computing Q_k(w)
    #####
    tabQLeg = struct_tabLeg.tabQLeg # Name of the container for Q_k(w)
    ######
    # Initial values of Q_k(w)
    val_0_Q = log(1.0-omg) - log(-1.0-omg) # Initial value for k = 0
    val_1_Q = 2.0 + omg*val_0_Q           # Initial value for k = 1
    #####
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg) # Computing the tabQLeg
    ##########
    # Computing D_k(w)
    #####
    for k=0:(K_BL-1) # Loop over the Legendre indices
        tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
    end
end
##################################################
# Function to fill in all the Legendre arrays
# for an NEUTRAL mode, i.e. Im[w] = 0.0
##################################################
function tabLeg!_NEUTRAL(omg::Complex{Float64},
                         struct_tabLeg::structLegendre_type)
    ##########
    tabDLeg = struct_tabLeg.tabDLeg # Name of the container for D_k(w)
    ##########
    romg = real(omg) # Keeping only the real part of omg
    ##########
    # Computing Q_k(w)
    #####
    tabQLeg = struct_tabLeg.tabQLeg # Name of the container for Q_k(w)
    ######
    # Initial values of Q_k(w)
    # The integral is computed in the Cauchy sense, hence the absolute values.
    val_0_Q = log(abs(1.0-romg)) - log(abs(-1.0-romg)) + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex. ATTENTION, not to forget the `abs'.
    val_1_Q = 2.0 + romg*val_0_Q                       + 0.0*im # Initial value for k = 1. ATTENTION, we make sure that the quantity is seen as complex.
    #####
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg) # Computing the tabQLeg
    ##########
    # Computing P_k(w) if needed
    #####
    val_H = Heaviside(romg) # Computing the value of H[-1,1,Re[w]]. ATTENTION, the argument must be real.
    #####
    if (val_H == 0.0) # We are outside of the interval [-1,1]: no need to compute P_k(w)
        ##########
        # Computing D_k(w)
        ##########
        for k=0:(K_BL-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays.
        end
        #####
    else # We are within the interval [-1,1]: we need to compute P_k(w)
        #####
        tabPLeg = struct_tabLeg.tabPLeg # Name of the container for P_k(w)
        #####
        # Initial values of P_k(w)
        val_0_P = 1.0  + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex.
        val_1_P = romg + 0.0*im # Initial value for k = 1. ATTENTION, we make sure that the quantity is seen as complex.
        #####
        tabPLeg!(omg,val_0_P,val_1_P,tabPLeg) # Computing the tabPLeg
        ##########
        # Computing D_k(w)
        #####
        for k=0:(K_BL-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] + im*pi*val_H*tabPLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
        end
        #####
    end
end
##################################################
# Function to fill in all the Legendre arrays
# for a DAMPED mode, i.e. Im[w] < 0.0
##################################################
function tabLeg!_DAMPED(omg::Complex{Float64},
                        struct_tabLeg::structLegendre_type)
    ##########
    tabDLeg = struct_tabLeg.tabDLeg # Name of the container for D_k(w)
    ##########
    romg = real(omg) # Keeping only the real part of omg
    ##########
    # Computing Q_k(w)
    #####
    tabQLeg = struct_tabLeg.tabQLeg # Name of the container for Q_k(w)
    ######
    # Initial values of Q_k(w)
    val_0_Q = log(1.0-omg) - log(-1.0-omg) # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex. ATTENTION, not to forget the `abs'.
    val_1_Q = 2.0 + omg*val_0_Q            # Initial value for k = 1. ATTENTION, we make sure that the quantity is seen as complex.
    #####
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg) # Computing the tabQLeg
    ##########
    # Computing P_k(w), if needed
    #####
    val_H = Heaviside(romg) # Computing the value of H[-1,1,Re[w]]. ATTENTION, the argument must be real.
    #####
    if (val_H == 0.0) # We are outside of the interval [-1,1], so that we do not need to compute P_k(w)
        ##########
        # Computing D_k(w)
        #####
        for k=0:(K_BL-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays.
        end
        #####
    else # We are within the interval [-1,1], so that we need to compute P_k(w)
        #####
        tabPLeg = struct_tabLeg.tabPLeg # Name of the container for P_k(w)
        #####
        # Initial values of P_k(w)
        val_0_P = 1.0  + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex.
        val_1_P = omg           # Initial value for k = 1.
        #####
        tabPLeg!(omg,val_0_P,val_1_P,tabPLeg) # Computing the tabPLeg
        ##########
        # Computing the D_k(w)
        #####
        for k=0:(K_BL-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] + 2.0*im*pi*val_H*tabPLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
        end
        #####
    end
end
##################################################
# Choosing the function tabLeg! depending
# on the type of linear response considered
##################################################
if (ALEPH == "unstable") # Searching for unstable modes, i.e. Im[w] > 0
    const tabLeg! = tabLeg!_UNSTABLE # ATTENTION, use `const' to avoid allocations
elseif (ALEPH == "neutral") # Searching for neutral modes, i.e. Im[w] = 0
    const tabLeg! = tabLeg!_NEUTRAL # ATTENTION, use `const' to avoid allocations
elseif (ALEPH == "damped") # Searching for damped modes, i.e. Im[w] < 0
    const tabLeg! = tabLeg!_DAMPED # ATTENTION, use `const' to avoid allocations
end
##################################################
# Function that pre-computes the Hilbert-transformed
# Legendre functions for a given complex frequency
# Q_k(w) = INT[P_k(u)/(u-w),{u,-1,1}]
# ATTENTION, this is == -2 q_k(w)
# with q_k(w) the Legendre functions of the second kind
# Arguments are:
# + omg: COMPLEX frequency.      ATTENTION, has to be complex.
# + val_0: Initial value in k=0. ATTENTION, has to be complex.
# + val_1: Initial value in k=1. ATTENTION, has to be complex.
# + tabQLeg: Container where to store the results
# There are two possible algorithms:
# + Close to the real line [-1,1], we use an upward recurrence
# + Far-away from the real line[-1,1], we use a downward recurrence
# The transition from the two regimes follows from the thesis
# Stable Implementation of Three-Term Recurrence Relations
# Pascal Frederik Heiter, June, 2010
# https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.070/funken/bachelorarbeiten/bachelorthesis_pfh.pdf
##################################################
function tabQLeg!(omg::Complex{Float64},
                  val_0::Complex{Float64},
                  val_1::Complex{Float64},
                  tabQLeg::Array{Complex{Float64},1})
    #####
    logINVTOL = log(10^(14)) # Logarithm of the inverse of the tolerance of the algorithm
    #####
    u, v = real(omg), imag(omg) # Real/Imag part of the frequency
    u2, v2 = u^(2), v^(2) # Squared quantities
    #####
    b = min(1.0,4.5/((K_BL+1.0)^(1.17))) # Value of ellise parameter
    b2 = b^(2) # Squared quantity
    a2 = 1.0+b2 # Second parameter of the ellipse
    #####
    # Testing whether or not we are close to the axis [-1,1]
    if ((u2/a2) + (v2/b2) <= 1.0) # Within the ellipse
        #####
        tabLeg_UP!(omg,val_0,val_1,tabQLeg) # If we are sufficiently close to [-1,1], we use the forward recurrence
    else # Outside of the ellipse
        #####
        # Determining the size of the warm-up region
        z2 = abs(u) + abs(v)*im
        ctmp0 = sqrt(z2*z2 - 1.0)
        tmp00 = 2.0*log(abs(z2 + ctmp0))
        #####
        K_c = convert(Int64,K_BL + ceil(logINVTOL/tmp00)) # Index at which the backward warm-up is started
        #####
        tabLeg_BACK!(omg,val_0,K_c,tabQLeg) # If we are sufficiently far away from [-1,1], we use the backward recurrence
    end
end
##################################################
# Function that pre-computes the Legendre functions, P_k(w),
# for a given complex frequency
# Arguments are:
# + omg: COMPLEX frequency. ATTENTION, has to be complex.
# + val_0: Initial value in k=0. ATTENTION, has to be complex.
# + val_1: Initial value in k=1. ATTENTION, has to be complex.
# + tabPLeg: Container where to store the results
# For the Legendre functions, we always use the forward recurrence
##################################################
function tabPLeg!(omg::Complex{Float64},
                  val_0::Complex{Float64},
                  val_1::Complex{Float64},
                  tabPLeg::Array{Complex{Float64},1})
    #####
    tabLeg_UP!(omg,val_0,val_1,tabPLeg) # For P_k(w), we always use the forward recurrence.
end
##################################################
# Function to compute Legendre functions
# with an UPWARD recurrence
# This function is used for both P_k(w) and Q_k(w)
# as they satisfy the same two-term recurrence
# + It is always used for P_k(w)
# + It is used for Q_k(w) only sufficiently close to [-1,1]
# Arguments are:
# + omg: COMPLEX frequency. ATTENTION, has to be complex.
# + val_0: Initial value for k=0. ATTENTION, has to be complex.
# + val_1: Initial value for k=1. ATTENTION, has to be complex.
# + tabLeg: Container where to store the result
##################################################
function tabLeg_UP!(omg::Complex{Float64},
                    val_0::Complex{Float64},
                    val_1::Complex{Float64},
                    tabLeg::Array{Complex{Float64},1})
    #####
    tabLeg[0+1] = val_0 # Filling in the value of D_0(omg). ATTENTION, to the shift in the array index
    tabLeg[1+1] = val_1 # Filling in the value of D_1(omg). ATTENTION, to the shift in the array index
    #####
    v0, v1 = val_0, val_1 # Initialising the temporary variables used in the recurrence
    #####
    for k=2:(K_BL-1) # Loop over the index 2 <= k
        v = ((2.0*k-1.0)*omg*v1 - (k-1.0)*v0)/(k) # Bonnet's recc. relation to get D_k(omg)
        #####
        tabLeg[k+1] = v # Filling in the value of D_k(omg). ATTENTION, to the shift in the array index
        #####
        v0, v1 = v1, v # Shifting the temporary variables
    end
end
##################################################
# Function to compute Legendre functions
# with a BACKWARD recurrence
# This function is only used for Q_k(w)
# sufficiently far away from [-1,1].
# Arguments are:
# + omg: COMPLEX frequency. ATTENTION, has to be complex.
# + val_0: Initial value in k=0. ATTENTION, has to complex.
# + K_c: Value of k at which the warm-up starts
# + tabLeg: Container where to store the results
##################################################
function tabLeg_BACK!(omg::Complex{Float64},
                      val_0::Complex{Float64},
                      K_c::Int64,
                      tabLeg::Array{Complex{Float64},1})
    #####
    v0 = 1.0 # Initialisation of D_(K_c+1)[omg]
    v1 = 0.0 # Initialisation of D_(K_c+2)[omg]
    #####
    # Warm-up phase
    for k=K_c:-1:K_BL # ATTENTION, the step is `-1'.
        v = ((2.0*k+3.0)*omg*v0 - (k+2.0)*v1)/(k+1.0) # Computing D_(k)[omg]
        v0, v1 = v, v0 # Shifting the temporary variables
    end
    ####
    # Starting to store the values of D_(k)[omg]
    for k=(K_BL-1):-1:0 # ATTENTION, the step is `-1'.
        v = ((2.0*k+3.0)*omg*v0 - (k+2.0)*v1)/(k+1.0) # Computing D_(k)[omg]
        #####
        tabLeg[k+1] = v # Filling in the value of D_(k)[omg]. ATTENTION, to the shift in the array index
        #####
        v0, v1 = v, v0 # Shifting the temporary variables
    end
    #####
    pref = val_0/tabLeg[0+1] # Prefactor by which all the values have to be rescaled, so that it is equal to val_0 for k=0. ATTENTION, to the shift of the array.
    #####
    for k=0:(K_BL-1) # Rescaling all the computed values
        tabLeg[k+1] *= pref # Performing the rescaling. ATTENTION, to the shift in the array index
    end
end



#################################################
# Computing the response matrix in omega thanks to the 
# previous functions
#################################################

#################################################
# Computing the part of response matrix associated
# to the n-th frequency in the FT
#################################################

function M_n_legendre!(n::Int64,omg::Complex{Float64},
                        structLeg,
                        Mmat::Array{Complex{Float64},2}) 

    tabompart = structLeg.tabOmegapart
    OmM, OmE = tabompart[3], tabompart[4]
    
    if (n == 0)
        return nothing
    elseif n > 0
        varpi =  omg / (n * OmE) - OmM / OmE # \varpi definition for n > 0
    else
        varpi =  - omg / (n * OmE) + OmM / OmE # \varpi definition for n < 0 
    end

    tabLeg!(varpi,structLeg) # Computing the values of D_k in varpi

    apqnk = structLeg.tabapqnk
    tabD = structLeg.tabDLeg

    for k = 1:K_BL 

        #= @timeit to "mat" =# begin
            for q = 1:NBASISELMT
                for p = 1:NBASISELMT
                    Mmat[p,q] += apqnk[p,q,n+(NMAXRESMAT+1),k] * tabD[k]
                end 
            end
        end
    end
end

##################################################
# Defining the response matrix method and initialization to Legendre
# ATTENTION NOT TO FORGET `const'
# OTHERWISE THERE ARE UNNECESSARY ALLOCATIONS
##################################################
const M_n! = M_n_legendre!
const BL_fullinit! = Legendre_fullinit!
