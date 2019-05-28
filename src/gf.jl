"""
Galois fields

Currently only supports p^n finite fields over the integers 0:p-1 where n=1 and p is prime. Does not support any fancy divisors.
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
    else #Multiplicative inverse does not exist
        throw(DivideError)
    end
end

#Promote any binary operation with some other number into the field
promote_rule(::Type{S}, x::Type{GF{T,q}}) where {q, S, T} = promote_rule(GF{S,q}, x)
#Use widest internal representation
promote_rule(::Type{GF{S,q}}, ::Type{GF{T,q}}) where {q, S, T} = GF{promote_type(S, T), q}

Random.rand(rng::AbstractRNG, ::Random.SamplerType{GF{T,q}}) where {T,q} =
    GF(T, q, rand(rng, 0:q-1))

end #module