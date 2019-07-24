# SecureComputation.jl
Secure multiparty computation in Julia

Implements Shamir's secret sharing scheme to distribute a number, and the BGW88 protocol for calculating + and * over the resulting numeric representation.

## Example

Example from [Wikipedia](https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing#Solution)

```jl
using SecureComputation

F = GF{Int,1613} #Galois field (finite field, integers mod 1613, perform compute in native Int)

secret = F(1234)
x = shamir(6, secret, F[166, 94]) # == DistributedShares(1:6, F[1494, 329, 965, 176, 1188, 775])

unshamir(x) # == secret

#With 2 random coefficients, you only need 2+1 = 3 pieces to reconstruct the secret
unshamir(x[1:3]) # == secret
```

A more realistic example using the Intel hardware CSPRNG

```
F = GF{BigInt,1000000004191}
secret = F(1234567890)
x = shamir(RdRand(), 100, 39, secret)
unshamir(x) # == secret

#If I have just 39 pieces, I can also reconstruct the number
unshamir(x[1:39])  # == secret
```

An example showing a simple matrix computation

```
using LinearAlgebra, SecureComputation
F = GF{BigInt,1000000004191}
M = rand(RdRand(), F, 3, 3)
MM = shamir.(100, M)
Z=lu(MM)
istril(unshamir.(Z.L)) #true
```

