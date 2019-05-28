############################################################################
# Shamir's secret sharing scheme
# An example of a linear secret sharing scheme
##############################################################################

struct DistributedShares{SD<:SecretDomain, U<:AbstractVector{Int}, V<:AbstractVector{SD}}
    idxs::U
    vals::V
end

"""
Evaluate a polynomial
c0 + x*cs[1] + x^2*cs[2] + ...
using Horner's rule
"""
function evalpoly(x, c0, cs)
    t = length(cs)
    if t > 0
        f = cs[t]
        for i=t-1:-1:1
            f = f * x + cs[i]
        end
        f = f * x + c0
    else
        f = c0
    end
end

"""
Evaluate the value of the Lagrange interpolating polynomial at 0
"""
function lagrangeinterp0(nodes, vals)
    @assert size(nodes) == size(vals)
    k = size(vals, 1)
    SD = eltype(vals)
    L = zero(SD)
    for j=1:k
        xj = nodes[j]
        c = vals[j]
        for m = 1:k
            if m != j
                xm = SD(nodes[m])
                c *= xm/(xm-xj)
            end
        end
        L += c
    end
    L
end
"""
    shamir([rng], t, [deg], s, [rs]) -> DistributedShare

Shamir's secret sharing scheme

Take a secret s embedded in some SecretDomain and produce t shares

The degree of shring defaults to the maximal degree which permits multiplication

Optionally specify the random number generator and
other coefficients of the polynomial used for embedding; otherwise choose at random

"""
function shamir(t::Int, s::SD, rs) where SD <: SecretDomain
    vals = zeros(SD, t)
    for j=1:t
        f = evalpoly(j, s, rs)
        vals[j] = f
    end
    DistributedShares(1:t, vals)
end

#Default degree of sharing is the maximal degree which permits multiplication
shamir(t::Int, secret::SD) where SD <: SecretDomain =
    shamir(GLOBAL_RNG, t::Int, secret::SD)
shamir(rng::AbstractRNG, t::Int, secret::SD) where SD <: SecretDomain =
    shamir(rng, t, (t-1)÷2, secret)

shamir(t::Int, n::Int, secret::SD) where SD <: SecretDomain =
    shamir(GLOBAL_RNG, t::Int, n::Int, secret::SD)
shamir(rng::AbstractRNG, t::Int, n::Int, secret::SD) where SD <: SecretDomain =
    shamir(t, secret, rand(rng, SD, n-1))

unshamir(x::DistributedShares) = lagrangeinterp0(x.idxs, x.vals)

