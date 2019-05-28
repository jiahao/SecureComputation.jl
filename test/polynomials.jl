#Test
using Polynomials, Test
let
    SD = GF{982451653}
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

    @test SecureComputation.reducedegree(an .* bn, SD.(1:n), t) == cn
end
