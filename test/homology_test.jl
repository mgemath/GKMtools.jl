@testset "second homology and curve classes" begin
  P2 = projective_space(GKMGraph, 2)
  @test !ismutable(P2)
  @test !ismutable(P2.H2)
  @test fieldtype(typeof(P2), :H2) === GKM_H2
  H2 = GKM_second_homology(P2)
  @test H2 === GKM_second_homology(P2)
  @test rank(H2.H2) == 1
  beta = curve_class(P2, first(GKMtools.edges(P2)))
  @test chern_number(P2, beta) == 3
  @test is_effective(P2, beta)
  @test !is_effective(P2, -beta)

  # Generalized flags use their root data instead of the generic, potentially
  # enormous cycle-relation matrix. B2 also checks the coroot normalization.
  B2 = generalized_gkm_flag(root_system(:B, 2), [2])
  flag_H2 = GKM_second_homology(B2)
  @test rank(flag_H2.H2) == 1
  flag_classes = unique(curve_class(B2, e)[1] for e in GKMtools.edges(B2))
  @test sort(flag_classes) == [1, 2]
  @test all(
    e -> chern_number(B2, curve_class(B2, e)) == 3curve_class(B2, e)[1],
    GKMtools.edges(B2),
  )

  fan = StackyFan([2 0; 0 3; -1 -1], [[1, 2], [2, 3], [3, 1]])
  X = gkm_graph_of_orbifold_toric(fan)
  @test !ismutable(X)
  @test !ismutable(X.H2)
  @test fieldtype(typeof(X), :H2) === GKM_H2
  orbifold_H2 = GKM_second_homology(X)
  @test orbifold_H2 === GKM_second_homology(X)
  @test rank(orbifold_H2.H2) == 1
  classes = [curve_class(X, e)[1] for e in GKMtools.edges(X)]
  @test sort(abs.(classes)) == [1, 1, 2]
  @test all(e -> chern_number(X, curve_class(X, e)) == 3abs(curve_class(X, e)[1]), GKMtools.edges(X))
end
