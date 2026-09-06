@testset "Oscar small quantum cohomology ring" begin
  P1 = projective_space(GKMGraph, 1)
  beta = curve_class(P1, Oscar.Edge(1, 2))
  QH, classes, q = small_quantum_cohomology_ring(P1; degrees=[beta])
  h = only(classes[2])
  q1 = only(q)

  @test is_graded(QH)
  @test length.(classes) == [1, 1]
  @test h^2 == q1
  @test degree(h)[1] == 2
  @test degree(q1)[1] == 4

  QHcodim, classes_codim, q_codim = small_quantum_cohomology_ring(
    P1; degrees=[beta], degree_convention=:codimension,
  )
  @test degree(only(classes_codim[2]))[1] == 1
  @test degree(only(q_codim))[1] == 2
  @test is_graded(QHcodim)

  @test_throws ArgumentError small_quantum_cohomology_ring(
    P1; degrees=[beta], degree_convention=:unknown,
  )
end

@testset "Threaded unmarked integration" begin
  P1 = projective_space(GKMGraph, 1)
  beta = 2 * curve_class(P1, Oscar.Edge(1, 2))
  products = [
    [point_class(P1, 1)],
    [first_chern_class(P1), point_class(P1, 2)],
  ]

  sequential = gromov_witten_nomarks(P1, beta, products; show_bar=false)
  threaded = gromov_witten_nomarks(
    P1, beta, products; show_bar=false, threaded=true,
  )
  @test threaded == sequential
end

@testset "Optimized unmarked insertion preprocessing" begin
  P2 = projective_space(GKMGraph, 2)
  beta = curve_class(P2, Oscar.Edge(1, 2))
  point = point_class(P2, 1)

  # Repeated point classes reuse the cached tangent Euler class.
  @test point_class(P2, 1) == point
  @test isassigned(GKMtools.get_cohomology(P2).euler_classes, 1)

  @test gromov_witten_nomarks(
    P2, beta, [point, point, point]; show_bar=false,
  ) == 2gens(GKMtools.equivariant_coefficient_ring(P2))[1] -
       gens(GKMtools.equivariant_coefficient_ring(P2))[2] -
       gens(GKMtools.equivariant_coefficient_ring(P2))[3]
  @test iszero(gromov_witten_nomarks(
    P2, beta, [point, point, point]; show_bar=false, fast_mode=true,
  ))
end

@testset "Rational connection coefficients in edge factors" begin
  u, w = QQ(2), QQ(7)
  @test GKMtools._b(u, w, QQ(2)) == 1 // (w * (w - u) * (w - 2u))
  @test GKMtools._b(u, w, QQ(-3)) == (w + u) * (w + 2u)
  @test_throws ErrorException GKMtools._b(u, w, QQ(1) // 2)
end

@testset "Small quantum cohomology of the A2 flag variety" begin
  G = generalized_gkm_flag(root_system(:A, 2))
  beta = curve_class(G, Oscar.Edge(1, 2))

  QH, classes, q = small_quantum_cohomology_ring(G; degrees=[beta])

  @test is_graded(QH)
  @test length.(classes) == [1, 2, 2, 1]
  @test length(q) == 2

  sigma_s1, sigma_s2 = classes[2]
  sigma_s2s1, sigma_s1s2 = classes[3]
  @test sigma_s1^2 == sigma_s2s1 + q[1]
  @test sigma_s2^2 == sigma_s1s2

  basis = collect(Iterators.flatten(schubert_basis(G)))
  c1 = first_chern_class(G)
  symmetric_indices = [(j, k) for j in eachindex(basis) for k in j:length(basis)]
  marked = [
    ev(1, c1) * ev(2, basis[j]) * ev(3, basis[k])
    for (j, k) in symmetric_indices
  ]
  unmarked = [
    GKMClass[c1, basis[j], basis[k]]
    for (j, k) in symmetric_indices
  ]
  @test gromov_witten(
    G, beta, 3, marked; show_bar=false, fast_mode=true,
  ) == gromov_witten_nomarks(
    G, beta, unmarked; show_bar=false, fast_mode=true,
  )
end
