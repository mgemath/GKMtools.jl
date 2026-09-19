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

@testset "vector_bundle_O on products" begin
  for enlarged in (false, true)
    E = vector_bundle_O([1, 2], [[2, 3], [-1, 0]]; enlarge_torus=enlarged)
    G = baseof(E)
    @test num_vertices(G) == 6
    @test valency(G) == 3
    @test rank(E) == 2
    @test rank_torus(G) == 5 + (enlarged ? 2 : 0)
    b = gens(lattice(G))
    @test E.weights[1, :] == (enlarged ? b[6:7] : [b[1] + b[3], b[1] + b[3]])
    # Check degrees on every invariant curve, including away from the origin.
    for v2 in 1:3
      u, v = 1 + 2*(v2-1), 2 + 2*(v2-1)
      @test E.weights[u, 1] - E.weights[v, 1] == 2 * weight(G, Edge(u, v))
      @test E.weights[u, 2] - E.weights[v, 2] == -weight(G, Edge(u, v))
    end
    for v1 in 1:2, a in 1:3, b in a+1:3
      u, v = v1 + 2*(a-1), v1 + 2*(b-1)
      @test E.weights[u, 1] - E.weights[v, 1] == 3 * weight(G, Edge(u, v))
      @test E.weights[u, 2] == E.weights[v, 2]
    end
    single = vector_bundle_O([2], [[-3], [0], [5]]; enlarge_torus=enlarged)
    scalar = vector_bundle_O(2, [-3, 0, 5]; enlarge_torus=enlarged)
    @test [[w[i] for i in 1:rank_torus(baseof(single))] for w in single.weights] ==
      [[w[i] for i in 1:rank_torus(baseof(scalar))] for w in scalar.weights]
    triple = vector_bundle_O([1, 2, 1], [[2, 0, -3]]; enlarge_torus=enlarged)
    H = baseof(triple)
    @test num_vertices(H) == 12
    @test triple.weights[6, 1] - triple.weights[12, 1] == -3 * weight(H, Edge(6, 12))
  end
  @test_throws ArgumentError vector_bundle_O(Int[], [[1]])
  @test_throws ArgumentError vector_bundle_O([1, 0], [[1, 2]])
  @test_throws ArgumentError vector_bundle_O([1, 2], Vector{Int}[])
  @test_throws ArgumentError vector_bundle_O([1, 2], [[1]])
  @test_throws ArgumentError vector_bundle_O([1, 2], [[1, 2]]; small_torus=true)
end
