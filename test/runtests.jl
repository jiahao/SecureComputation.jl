using SecureComputation, Test

@testset "GF Test" begin include("gf.jl") end
@testset "Polynomials Test" begin include("polynomials.jl") end
@testset "Shamir Test" begin include("shamir.jl") end
