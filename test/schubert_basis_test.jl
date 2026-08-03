@testset "Schubert basis via Billey-style construction" begin
  R = root_system(:A, 2)
  G = generalized_gkm_flag(R)

  basis = schubert_basis(G)
  @test length.(basis) == [1, 2, 2, 1]
  @test sum(length, basis) == GKMtools.num_vertices(G)
  @test all(b -> parent(b) == GKMtools.get_cohomology(G), Iterators.flatten(basis))
  @test length.(billey_schubert_basis(G)) == [1, 2, 2, 1]
  @test schubert_basis(G, "s1*s2") isa Any

  polynomial_basis = schubert_basis(G; representation=:polynomial)
  @test length.(polynomial_basis) == [1, 2, 2, 1]
  @test all(c -> c isa GKMClass, Iterators.flatten(polynomial_basis))
  @test localize.(Iterators.flatten(polynomial_basis)) == collect(Iterators.flatten(basis))
  @test_throws ArgumentError schubert_basis(G; representation=:unknown)
end
