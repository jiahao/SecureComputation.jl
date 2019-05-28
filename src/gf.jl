"""
Galois fields
"""
module GaloisFields

export SecretDomain, FiniteField, GF
using Random
import Base: +, -, *, /, inv, promote_rule, rand

abstract type SecretDomain <: Number end
abstract type FiniteField <: SecretDomain end

struct GF{q} <: FiniteField
    val :: Int
end

GF(q, a) = GF{q}(a % q)

+(a::GF{q}, b::GF{q}) where q = GF{q}(mod(a.val+b.val, q))
-(a::GF{q}, b::GF{q}) where q = GF{q}(mod(a.val-b.val, q))
-(a::GF{q}) where q = GF{q}(mod(-a.val, q))
*(a::GF{q}, b::GF{q}) where q = GF{q}(mod(a.val*b.val, q))
/(a::GF{q}, b::GF{q}) where q = a * inv(b)

function inv(a::GF{q}) where q
    _, ainv, _ = gcdx(a.val, q)
    return GF{q}(mod(ainv, q))
end

#Promote any binary operation with some integer into the field
promote_rule(::Type{GF{q}}, ::Type{T}) where q where T<:Integer = GF{q}
promote_rule(::Type{T}, ::Type{GF{q}}) where q where T<:Integer = GF{q}

rand(rng::AbstractRNG, ::Random.SamplerType{GF{q}}) where q =
    GF{q}(rand(rng, 0:q-1))
end #module

using .GaloisFields, Test
let
    F = GF{1613}

    @test F(1) + F(1) == F(2)
    @test F(1) - F(1) == F(0)
    @test F(1) * F(1) == F(1)
    @test F(1) / F(1) == F(1)
    @test inv(F(1)) == F(1)

    @test F(1) + 1 == F(2)
    @test 1 + F(1) == F(2)

    @test typeof(rand(F, 5)) == Vector{F}
end

