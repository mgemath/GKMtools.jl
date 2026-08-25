@testset "vector_bundle_O torus enlargement" begin
  degrees = [-3, 0, 5]

  current = vector_bundle_O(3, degrees)
  @test rank_torus(baseof(current)) == 4
  @test rank(current) == 3

  former = vector_bundle_O(3, degrees; enlarge_torus=true)
  @test rank_torus(baseof(former)) == 7
  @test rank(former) == 3
  @test former.weights[1, :] == gens(lattice(baseof(former)))[5:7]

  line = vector_bundle_O(2, [1]; enlarge_torus=true)
  P2 = baseof(line)
  beta = curve_class(P2, "1", "2")
  M = twisted_c1_matrix(line, beta)
  t = gens(parent(numerator(M[1, 1])))
  expected = 2 * t[4] // ((t[1] - t[2]) * (t[1] - t[3]))
  @test all(==(expected), M[1, :])
  expected_second =
    (2*t[1] - 2*t[2] - 2*t[4]) // ((t[1] - t[2]) * (t[2] - t[3]))
  expected_third =
    (-2*t[1] + 2*t[3] + 2*t[4]) // ((t[1] - t[3]) * (t[2] - t[3]))
  @test all(==(expected_second), M[2, :])
  @test all(==(expected_third), M[3, :])
end
