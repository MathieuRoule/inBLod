const binjl = true  # Flag of file reading
###################################################
# Bin center
#
# Bin k (1:NBIN) of a quantity a : 
#   [a_min + (k-1)*da, a_min + k*da [
# Then, the bin centers are
#       a_min + (k-0.5)*da
###################################################
function bincenter(k::Int64,a_min::Float64,da::Float64)
    return a_min + (k - 0.5) * da
end

###################################################
# Adding each particle to the bins in which it is 
#
# Bin k (1:NBIN) of a quantity a : 
#   [a_min + (k-1)*da, a_min + k*da [
# Then, a in bin k means
#       a_min + (k-1)*da =< a < a_min + k*da
# <=>   k-1 =< (a-a_min)/da < k
# <=>   k = floor( (a-a_min)/da ) + 1
###################################################
function whichbin(a::Float64,a_min::Float64,da::Float64)
    return floor(Int64,(a-a_min)/da) + 1 # Could possibly be higher than NBINS --> to avoid in counting after
end
