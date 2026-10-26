@testset "GKM graphs-standard constructions" begin
  @testset "Flag variety axial weights and tautological bundles" begin
    G = flag_variety(GKMGraph, [1, 1])
    edge = only(edges(G))
    basis = gens(lattice(G))

    @test weight(G, edge) == basis[1] - basis[2]

    S, _ = tautological_and_univ_bd(GKMGraph, 1, 2)
    dual_S = dual(S)
    @test GKMtools.coefficients(S.connection)[edge] == [-1]
    @test GKMtools.coefficients(dual_S.connection)[edge] == [1]

    L = line_bundle_O(1, 1)
    @test GKMtools.coefficients(L.connection)[edge] == [1]
  end
end