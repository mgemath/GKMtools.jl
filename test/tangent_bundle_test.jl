@testset "Tangent bundle fibre scaling" begin
  G = projective_space(GKMGraph, 2)
  E = tangent_bundle(G)
  @test rank(E) == 2
  @test rank(E.M) == 4
  @test GKMtools.transport(GKMtools.connection(E)) == GKMtools.transport(GKMtools.connection(G))
  expected = [
    (1, -1, 0, 1) (1, 0, -1, 1);
    (-1, 1, 0, 1) (0, 1, -1, 1);
    (-1, 0, 1, 1) (0, -1, 1, 1)
  ]
  for v in 1:3, i in 1:2
    @test fiber_weight(E, v, i) == E.M(collect(expected[v, i]))
  end

  unscaled = tangent_bundle(G; big_torus=false)
  @test unscaled.M === lattice(G)
  @test all(fiber_weight(unscaled, v, i) == weight(G, v, i) for v in 1:3, i in 1:2)
  @test integrate(chern_class(G, 2)) == 3
  @test integrate(first_chern_class(G)^2) == 9
  @test chern_classes(G)[3] == chern_class(G, 2)
  @test total_chern_class(G) == sum(chern_classes(G))

  fan = StackyFan([1 0; 0 1; -1 -1], [[1, 2], [2, 3], [3, 1]]; generic_stabilizer=[2])
  orbifold = gkm_graph_of_orbifold_toric(fan)
  scaled = tangent_bundle(orbifold)
  original = tangent_bundle(orbifold; big_torus=false)
  @test rank(scaled.M) == rank(lattice(orbifold)) + 1
  for v in 1:num_vertices(orbifold), i in 1:rank(scaled)
    @test fiber_weight(scaled, v, i) == scaled.GMtoM(fiber_weight(original, v, i)) + gens(scaled.M)[end]
    @test fiber_representation(scaled, v) == fiber_representation(original, v)
  end
end
