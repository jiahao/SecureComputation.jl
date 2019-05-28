module SecureComputation

export shamir, unshamir, GF, RdRand, RdSeed

using Random
GLOBAL_RNG = Random.GLOBAL_RNG

include("cpurng.jl")
using .CPURNGs

include("gf.jl")
include("shamir.jl")

import Base: +, -, *, /, \, inv, eltype, zero, oneunit, promote_type, getindex

#Array operations on DistributedShares

getindex(x::DistributedShares, inds) = DistributedShares(inds, x.vals[inds])

#Arihmetic on DistributedShares

zero(x::DistributedShares) = DistributedShares(x.idxs, x.vals*0)
zero(x::Type{DistributedShares{T,S,R}}) where {T,S,R} = DistributedShares(1:0, T[])

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
+(a::T, b::T, ::Val{true}) where T <: DistributedShares = DistributedShares(a.idxs, a.vals + b.vals)

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
    DistributedShares(idxs, vals)
end

##############################################################################
# Subtraction
##############################################################################

-(a::T) where T <: DistributedShares = DistributedShares(a.idxs, -a.vals)
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
function randomize!(x::DistributedShares, t::Int)
    vals = x.vals
    SD = eltype(x)
    n = length(vals)
    for i in x.idxs
        q = shamir(n, t, SD(0))
        vals += q.vals
    end
    DistributedShares(x.idxs, vals)
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

#Test
using Polynomials
let
    SD = GF{BigInt,982451653}
    #Number of interpolating points
    n = 11
    #Degree
    t = 4


    Pa = Poly(SD.(1:5))
    Pb = Poly(SD[6,7,8,9,1])
    #Truncate the product to degree t
    Pc = Poly((Pa * Pb).a[1:t+1])

    an = [polyval(Pa, i) for i=1:n]
    bn = [polyval(Pb, i) for i=1:n]
    cn = [polyval(Pc, i) for i=1:n]

    @test reducedegree(an .* bn, SD.(1:n), t) == cn
end

*(a::T, b::T) where T <: DistributedShares =
    *(a::T, b::T, Val(isalignedindexes(a, b)))

function *(a::T, b::T, ::Val{false}) where T <: DistributedShares
    S = eltype(T)
    idxs = alignindexes(a.idxs, b.idxs)
    n = length(idxs)
    #XXX HARD CODED ASSUMPTION
    # SHOULD PUT T IN TYPE PARAMETER
    # n = 2t + 1
    t = (n - 1) ÷ 2

    vals = zeros(S, n)
    for i in (a.idxs) ∩ (b.idxs)
        vals[i] = a.vals[i] * b.vals[i]
    end

    #Degree reduction
    R = reducedegree(vals, SD.(idxs), t)
    h = DistributedShares(idxs, R)

    #Randomization
    randomize!(h, t)
end

function *(a::T, b::T, ::Val{true}) where T <: DistributedShares
    n = length(a.vals)
    idxs = a.idxs

    #XXX HARD CODED ASSUMPTION
    # SHOULD PUT T IN TYPE PARAMETER
    # n = 2t + 1
    t = (n - 1) ÷ 2

    SD = eltype(T) #get SecretDomain

    #Degree reduction
    R = reducedegree(a.vals.*b.vals, SD.(idxs), t)
    h = DistributedShares(idxs, R)

    #Randomization
    randomize!(h, t)
end

##############################################################################
# Inversion
##############################################################################

inv(a::DistributedShares) = DistributedShares(a.idxs, inv.(a.vals))

##############################################################################
# Division
##############################################################################

/(a::T, b::T) where T <: DistributedShares = a * inv(b)
\(a::T, b::T) where T <: DistributedShares = inv(a) * b

##############################################################################
using .Test

let
    SD = GF{BigInt,982451653}
    secret1 = 12345
    secret2 = 54321

    n = 16
    t = 6

    x = shamir(n, t, SD(secret1))
    y = shamir(n, t, SD(secret2))

    @test unshamir(x+y) == SD(secret1+secret2)
    @test unshamir(x-x) == SD(0)
    @test unshamir(x+(-x)) == SD(0)
    @test unshamir(x*y) == SD(secret1*secret2)
    @test unshamir(x/x) == SD(1)
    @test unshamir(x\x) == SD(1)
end

let
    SD = GF{BigInt,982451653}
    n = 4
    t = 6
    x = rand(SD, n, n)
    y = rand(SD, n, n)
    z = x + y

    sx = shamir.(t, t, x)
    sy = shamir.(t, t, y)
    sz = sx + sy

    @test unshamir.(sz) == z
end

let
    SD = GF{BigInt,982451653}
    n = 1
    t = 3
    x = rand(SD, n, n)
    y = rand(SD, n, n)
    z = x * y

    sx = shamir.(2t+1, t, x)
    sy = shamir.(2t+1, t, y)
    sz = sx * sy

    @test unshamir.(sz) == z
end

end # module

