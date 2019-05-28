"""
Galois fields
"""
module GaloisFields

export SecretDomain, FiniteField, GF
using Random
import Base: ==, +, -, *, /, \, oneunit, inv, promote_rule, rand, convert

abstract type SecretDomain <: Number end
abstract type FiniteField <: SecretDomain end

struct GF{T,q} <: FiniteField
    val :: T
end

GF(q, a) = GF(typeof(a), q, a)
GF(T::Type, q, a) = GF{T,q}(mod(T(a), T(q)))

oneunit(::Type{GF{T,q}}) where {T, q} = GF{T,q}(1)
convert(::Type{GF{T,q}}, x::GF{S,q}) where {S,T,q} = GF{T,q}(x.val)

+(a::GF{T,q}, b::GF{T,q}) where {T,q} = GF{T,q}(mod(a.val+b.val, q))
-(a::GF{T,q}, b::GF{T,q}) where {T,q} = GF{T,q}(mod(a.val-b.val, q))
-(a::GF{T,q}            ) where {T,q} = GF{T,q}(mod(-a.val, q))
*(a::GF{T,q}, b::GF{T,q}) where {T,q} = GF{T,q}(mod(a.val*b.val, q))
/(a::GF{T,q}, b::GF{T,q}) where {T,q} = a * inv(b)
\(a::GF{T,q}, b::GF{T,q}) where {T,q} = inv(a) * b
==(a::GF{T,q}, b::GF{T,q}) where {T,q} = a.val == b.val
function inv(a::GF{T,q}) where {T,q}
    x, b, _ = gcdx(a.val, q)
    if x == 1
        return GF{T, q}(mod(b, q))
    else
        throw(DivideError)
    end
end

#XXX Missing promote_rule in base?
promote_rule(::Type{<:Integer}, ::Type{BigInt}) = BigInt

#Promote any binary operation with some integer into the field
promote_rule(x::Type{GF{T,q}}, ::Type{S}) where {q, S, T} = promote_rule(x, GF{S,q})
promote_rule(::Type{S}, x::Type{GF{T,q}}) where {q, S, T} = promote_rule(GF{S,q}, x)
promote_rule(::Type{GF{S,q}}, ::Type{GF{T,q}}) where {q, S, T} = GF{promote_rule(S, T), q}
promote_rule(::Type{GF{T,q}}, ::Type{GF{T,q}}) where {q, T} = GF{T, q}


Random.rand(rng::AbstractRNG, ::Random.SamplerType{GF{T,q}}) where {T,q} =
    GF(T, q, rand(rng, 0:q-1))

end #module

using .GaloisFields, Test
let
    F = GF{Int,1613}

    @test F(1) + F(1) == F(2)
    @test F(1) - F(1) == F(0)
    @test F(1) * F(1) == F(1)
    @test F(1) / F(1) == F(1)
    @test inv(F(1)) == F(1)

    @test F(1) + 1 == F(2)
    @test 1 + F(1) == F(2)


    @test typeof(rand(F)) == F
    @test typeof(rand(F, 5)) == Vector{F}
end

