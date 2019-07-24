module SecureComputation

export shamir, unshamir, GF, RdRand, RdSeed

using Random
GLOBAL_RNG = Random.GLOBAL_RNG

include("cpurng.jl")
using .CPURNGs

include("gf.jl")

using .GaloisFields

include("shamir.jl")

export GF, DistributedShares, shamir, unshamir

# Arithmetic on distributed shares
import Base: +, -, *, /, \, inv, eltype, zero, oneunit, promote_type, getindex

#Array operations on DistributedShares

getindex(x::DistributedShares{SD, N, D, U, V}, inds) where {SD, N, D, U, V} =
    DistributedShares(inds, D, x.vals[inds])

#Arithmetic on DistributedShares

#zero(x::DistributedShares{SD, N, D, U, V}) where {SD, N, D, U, V} =
#    DistributedShares(x.idxs, D, x.vals*0)
#zero(x::Type{DistributedShares{SD, N, D, U, V}}) where {SD, N, D, U, V} =
#    DistributedShares(1:0, D, T[])

eltype(::Type{S}) where S <: DistributedShares = S.parameters[1]

function promote_type(::Type{DistributedShares{Ta, Ua, Va}}, b::Type{DistributedShares{Tb, Ub, Vb}}) where {Ta, Tb, Ua, Ub, Va, Vb}
    T = promote_type(Ta, Tb)
    U = promote_type(Ua, Ub)
    V = promote_type(Va, Vb)
    DistributedShares{T,U,V}
end

"Are indexes of shares aligned?"
isalignedindexes(a::DistributedShares, b::DistributedShares) =
    a.idxs == b.idxs

##############################################################################
# ARITHMETIC ACCORDING TO THE BGW88 PROTOCOL
# BGW88 defines + and *
# -, \, /, inv are the "obvious" inverse operations
##############################################################################

##############################################################################
# Addition
##############################################################################

+(a::T, b::T) where T <: DistributedShares =
    +(a::T, b::T, Val(isalignedindexes(a, b)))

#Addition when indexes are aligned
+(a::T, b::T, ::Val{true}) where T <: DistributedShares = DistributedShares(a.idxs, T.parameters[3], a.vals + b.vals)

function alignindexes(a::T, b::T) where T<:UnitRange
    UnitRange(min(a.start, b.start), max(a.stop, b.stop))
end

#Addition when indexes are not aligned
function +(a::T, b::T, ::Val{false}) where T <: DistributedShares
    S = eltype(T)
    idxs = alignindexes(a.idxs, b.idxs)
    n = length(idxs)
    vals = zeros(S, n)
    for i in idxs
        if i in a.idxs
            vals[i] += a.vals[i]
        end
        if i in b.idxs
            vals[i] += b.vals[i]
        end
    end
    DistributedShares(idxs, T.parameters[3], vals)
end

##############################################################################
# Subtraction
##############################################################################

-(a::T) where T <: DistributedShares = DistributedShares(a.idxs, T.parameters[3], -a.vals)
-(a::T, b::T) where T <: DistributedShares = a + (-b)

##############################################################################
# Multiplication
# Source: BGW88
# Completeness Theorems for Non-Cryptographic Fault-Tolerant Distributed Computation
# Michael Ben-Or and Shafi Goldwasser and Avi Wigderson
# 1988
##############################################################################

struct Vandermonde{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end

import Base: size, \, getindex, Matrix

oneunit(::Type{GF{T}}) where T = GF{T}(1)

size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(V::Vandermonde) = length(V.c), length(V.c)

function Matrix(V::Vandermonde{T}) where T
    n=size(V, 1)
    M=Array{T}(undef, n, n)
    M[:,1] .= 1
    for j=2:n
        for i=1:n
	    M[i,j] = M[i,j-1]*V.c[i]
        end
    end
    M
end

function vandtype(T1::Type, T2::Type)
    # Figure out the return type of Vandermonde{T1} \ Vector{T2}
    T = promote_type(T1, T2)
    S = typeof(oneunit(T)/oneunit(T1))
    return S
end

function \(V::Vandermonde{T1}, y::AbstractVecOrMat{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, size(y))
    copyto!(x, y)
    dvand!(V.c, x)
    return x
end

### HACK TO GET MATMUL
getindex(M::Vandermonde, i::Int, j::Int) = (M.c[i])^(j-1)

"""
    dvand!(a, b) -> b
Solves system ``A*x = b`` in-place.
``A`` is Vandermonde matrix ``A_{ij} = a_i^{j-1}``.
Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function dvand!(alpha, B)
    n = length(alpha)
    if n != size(B,1)
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
    end
    nrhs = size(B,2)
    @inbounds begin
        for j=1:nrhs
            for k=1:n-1
                for i=n:-1:k+1
                    B[i,j] = (B[i,j]-B[i-1,j])/(alpha[i]-alpha[i-k])
                end
            end
            for k=n-1:-1:1
                for i=k:n-1
                    B[i,j] = B[i,j]-alpha[k]*B[i+1,j]
                end
            end
        end
    end
end

#Add n random polynomials of degree t
function randomize!(x::T, t::Int) where T<:DistributedShares
    vals = x.vals
    SD = eltype(x)
    n = length(vals)
    for i in x.idxs
        q = shamir(n, t, SD(0))
        vals += q.vals
    end
    DistributedShares(x.idxs, T.parameters[3], vals)
end

#NOTE: Notation is transposed relative to the BGW88 paper,
# which operates on row vectors
function reducedegree(S, nodes, t::Integer)
    Te = eltype(S)
    n = length(S)

    B = Vandermonde(nodes)
    H = B \ S
    K = H
    for i=t+2:n
        K[i] = 0
    end
    B * K
end

*(a::T, b::T) where T <: DistributedShares =
    *(a::T, b::T, Val(isalignedindexes(a, b)))

function *(a::T, b::T, ::Val{false}) where T <: DistributedShares
    S = eltype(T)
    idxs = alignindexes(a.idxs, b.idxs)
    n = length(idxs)
    t = T.parameters[3]

    vals = zeros(S, n)
    for i in (a.idxs) ∩ (b.idxs)
        vals[i] = a.vals[i] * b.vals[i]
    end

    #Degree reduction
    R = reducedegree(vals, SD.(idxs), t)
    h = DistributedShares(idxs, t, R)

    #Randomization
    randomize!(h, t)
end

function *(a::T, b::T, ::Val{true}) where T <: DistributedShares
    n = length(a.vals)
    idxs = a.idxs

    SD = eltype(T) #get SecretDomain
    t = T.parameters[3]

    #Degree reduction
    R = reducedegree(a.vals.*b.vals, SD.(idxs), t)
    h = DistributedShares(idxs, t, R)

    #Randomization
    randomize!(h, t)
end

##############################################################################
# Inversion
##############################################################################

inv(a::T) where T<:DistributedShares = DistributedShares(a.idxs, T.parameters[3], inv.(a.vals))

##############################################################################
# Division
##############################################################################

/(a::T, b::T) where T <: DistributedShares = a * inv(b)
\(a::T, b::T) where T <: DistributedShares = inv(a) * b

##############################################################################

end # module
