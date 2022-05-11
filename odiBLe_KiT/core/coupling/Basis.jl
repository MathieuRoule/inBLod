const L_per = L_p * Lambda # Basis periodization lengthscale
##################################################
# Structure used to compute the bare susceptibility
# coefficients using the basis
##################################################
struct structBasis_type
    tabpsicos ::Array{Float64,1} # Basis cosine elmts value on a specific point
    tabpsisin ::Array{Float64,1} # Basis sine elmts value on a specific point
    tabftpsicos ::Array{Float64,1} # Fourier transform of each cos basis elmt for k
    tabftpsisin ::Array{Float64,1} # Fourier transform of each sin basis elmt for k
    tabftpsicosp ::Array{Float64,1} # Fourier transform of each cos basis elmt for k'
    tabftpsisinp ::Array{Float64,1} # Fourier transform of each sin basis elmt for k'
end
##################################################
# Creates a Basis_type
##################################################
function structBasis_create()
    return structBasis_type(zeros(Float64,NBASISELMT),zeros(Float64,NBASISELMT), # cosine and sine
                            zeros(Float64,NBASISELMT),zeros(Float64,NBASISELMT), # Fourier transform for k
                            zeros(Float64,NBASISELMT),zeros(Float64,NBASISELMT)) # Fourier transform for k'
end


##################################################
# For the cos/sin basis (periodized potential)
# Function to calculate the value of the basis elements 
# at a specific point. We do all the vector of the basis at the 
# same time using recurrence relation .
# Using Numerical Recipes eq (5.4.6)
# for sin/cos ((2p+1)y) calculation using recurrence. 
##################################################
function basisvalue_cossin!(x::Float64,tabpcos::Array{Float64,1},tabpsin::Array{Float64,1})

    pref = 2.0 * sqrt(G * L_per) / pi # Same prefactor for all the basis elements
    y = pi * x / L_per # Basis propto sin or cos ( (2p+1)*y )

    tabpsin[1], tabpcos[1] = sincos(y) # First basis sine/cosine element
    a = 2.0*(sin(y))^(2) # Used for calculation 
    b = sin(2.0*y) # Used for calculation 

    # We deduce the value for 2*(p+1)+1 using the value in 2 and 2*p+1
    for p=2:NBASISELMT 
        @inbounds c = tabpcos[p-1]
        @inbounds s = tabpsin[p-1]

        @inbounds tabpcos[p] = c - (a*c + b*s) # cos((2*(p+1)+1)y) = cos((2p+1)y) - (a*cos((2p+1)y) + b*sin((2p+1)y))
        @inbounds tabpsin[p] = s - (a*s - b*c) # sin((p+1)y) = sin((2p+1)y) - (a*sin((2p+1)y) - b*cos((2p+1)y))
    end
    for p=1:NBASISELMT # Prefactor of the (co)sine in the basis elements
        @inbounds tabpcos[p] *= pref / (2.0*p-1.0) # 2p-1 because p = 1,...,N (first element of the basis p=0)
        @inbounds tabpsin[p] *= pref / (2.0*p-1.0) # idem
    end
end

#####
# Possible others basis elements but only CosSin basis at the moment
# tabcos (resp. tabsin) should be used as tabs for even (resp. odd) basis elements 
# if another basis is built.
#####
if BASIS == "CosSin"
    const basisvalue! = basisvalue_cossin!
end

##################################################
# Function to calculate the Fourier transformed basis 
# element. We do all the vector of the basis at the 
# same time !
# If k is even the fourier transform (FT) of sine basis elmt is equal to 0.
# If k is odd the FT of cosine basis elmt is equal to 0.
##################################################
function fourierbasisvalue!(k::Int64,ra::Float64,
                   tabpcos::Array{Float64,1},tabpsin::Array{Float64,1}, # Basis elmt value
                   tabftpcos::Array{Float64,1},tabftpsin::Array{Float64,1}) # Basis Fourier transform elmt value
    #####
    Delta_u = 2.0/(K_FT) # Considered step distance
    Omega = Omegafromra(ra) # Value of the frequency
    #####
    for p = 1:NBASISELMT
        @inbounds tabftpcos[p] = 0.0 # Initialisation of the integrals 
        @inbounds tabftpsin[p] = 0.0 # Initialisation of the integrals
    end
    ##### 
    # Warm-up for the first point
    #####
    u = -1.0 
    x = ra*fu(u) # Initial radius of the particle
    dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du 
    # Step 1
    k1 = 0.5*Delta_u*dtdu # Current velocity for theta
    # Step 2
    u += 0.25*Delta_u # Updating the time by half a timestep
    x = ra*fu(u) # Current location of the radius, r=r(u)
    dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du
    #####
    k2 = 0.5*Delta_u*dtdu # Current velocity for theta
    # Step 3
    k3 = k2 # Does not need to be updated
    ########################################
    # Step 4
    # First true point
    u += 0.25*Delta_u # Updating the time by half a timestep
    x = ra*fu(u) # Current location of the radius, r=r(u)
    dtdu = Omega * dthetadu_omega(ra,u) # Current value of dtheta1/du and dtheta2/du, always well-posed
    #####
    k4 = 0.5*Delta_u*dtdu # Current velocity for theta
    #####
    # Better guess for the value of theta(u=-1.0+Delta_u/2.0)
    # with a first RK4 warm-up integration for a small step
    theta = (k1 + 2.0*k2 + 2.0*k3 + k4)/(6.0) # Updating theta
    #####
    basisvalue!(x,tabpcos,tabpsin) # Value of the basis element in x
    fact = dtdu*cos(k*theta)
    if (mod(k,2) == 0) # k is even, FT of sine elmt is 0
        for p = 1:NBASISELMT
            @inbounds tabftpcos[p] += tabpcos[p]*fact # First point u_1 contribution to the integral
        end
    else # k is odd, FT of cosine elmt is 0
        for p = 1:NBASISELMT
            @inbounds tabftpsin[p] += tabpsin[p]*fact # First point u_1 contribution to the integral
        end
    end
    #####
    for istep=2:K_FT # Loop over the steps of the integration !! ATTENTION, we only make (K_FT-1) steps
        #####
        # Integrating theta for one timestep
        ########################################
        # Step 1
        k1 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Step 2
        u += 0.5*Delta_u # Updating the time by half a timestep
        x = ra*fu(u) # Current location of the radius, r=r(u)
        dtdu = Omega * dthetadu_omega(ra,u) # Value of dtheta/du
        #####
        k2 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Step 3
        k3 = k2 # Does not need to be updated
        ########################################
        # Step 4
        u += 0.5*Delta_u # Updating the time by half a timestep
        x = ra*fu(u) # Current location of the radius, r=r(u)
        dtdu = Omega * dthetadu_omega(ra,u) # Current value of dtheta1/du and dtheta2/du, always well-posed
        #####
        k4 = Delta_u*dtdu # Current velocity for theta
        ########################################
        # Update
        theta += (k1 + 2.0*k2 + 2.0*k3 + k4)/(6.0) # Updating theta
        #####
        # Putting the values of the radius and prefactors
        basisvalue!(x,tabpcos,tabpsin) # Value of the basis element in x
        fact = dtdu * cos(k*theta)
        if (mod(k,2) == 0) # k is even, FT of sine elmt is 0
            for p = 1:NBASISELMT 
                @inbounds tabftpcos[p] += tabpcos[p]*fact # Point u_istep contribution to the integral
            end
        else # k is odd, FT of cosine elmt is 0
            for p = 1:NBASISELMT 
                @inbounds tabftpsin[p] += tabpsin[p]*fact # Point u_istep contribution to the integral
            end   
        end  
    end
    pref = 2.0 / (pi * K_FT) # Normalisation of the sum approx of the integral
    if (mod(k,2) == 0) # k is even, FT of sine elmt is 0
        for p = 1:NBASISELMT 
            @inbounds tabftpcos[p] *= pref # Normalisation of the sum approx of the integral
        end 
    else # k is odd, FT of cosine elmt is 0
        for p = 1:NBASISELMT 
            @inbounds tabftpsin[p] *= pref # Normalisation of the sum approx of the integral
        end 
    end   
end
