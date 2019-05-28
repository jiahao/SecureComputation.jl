using SecureComputation.GaloisFields, Test
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

#Protocol Pseudorandom Secret-Sharing (PRSS)
let
    #Example from
    #https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing#Solution
    F = GF{1613}

    secret = F(1234)
    idxs = 1:6
    xs = F[1494, 329, 965, 176, 1188, 775]

    @test shamir(6, secret, F[166, 94]).vals == xs
    @test unshamir(DistributedShares(idxs, xs)) == secret

    @test unshamir(shamir(6, secret)) == secret
end
