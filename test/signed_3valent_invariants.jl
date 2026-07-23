@testset "Signed 3-valent GKM invariant systems" begin
  function manual_system(mu_values, c1_values, p1_values; b3=0, euler=0)
    r = length(c1_values)
    triples = GKMtools._symmetric_triples(r)
    @assert length(mu_values) == length(triples)
    @assert length(p1_values) == r
    packed = matrix(ZZ, 1, length(triples), ZZ.(mu_values))
    c1 = matrix(ZZ, r, 1, ZZ.(c1_values))
    p1 = matrix(ZZ, 1, r, ZZ.(p1_values))

    tensor_entry(i, j, k) = packed[1, findfirst(==(Tuple(sort([i, j, k]))), triples)]
    tensor_eval(x, y, z) = sum(
      tensor_entry(i, j, k) * x[i] * y[j] * z[k]
      for i in 1:r for j in 1:r for k in 1:r; init=ZZ(0)
    )
    l11 = matrix(ZZ, 1, r, [
      tensor_eval(c1_values, c1_values, [i == j ? 1 : 0 for i in 1:r]) for j in 1:r
    ])
    differences = [l11[1, i] - p1[1, i] for i in 1:r]
    @assert all(iseven, differences)
    c2 = matrix(ZZ, 1, r, [divexact(d, ZZ(2)) for d in differences])

    R, z = polynomial_ring(ZZ, ["z$i" for i in 1:r])
    cubic = sum(
      GKMtools._multiplicity(i, j, k) * packed[1, q] * z[i] * z[j] * z[k]
      for (q, (i, j, k)) in enumerate(triples); init=zero(R)
    )
    c1_cubed = tensor_eval(c1_values, c1_values, c1_values)
    p1_c1 = (p1 * c1)[1, 1]
    characteristic = (
      c1_cubed=c1_cubed,
      c1_c2=(c2 * c1)[1, 1],
      c3=ZZ(euler),
      p1_c1=p1_c1,
    )
    return GKMtools._GKMInvariantSystem(
      r,
      free_module(ZZ, r),
      zero_matrix(ZZ, 0, r),
      packed,
      triples,
      R,
      cubic,
      c1,
      change_base_ring(GF(2), c1),
      p1,
      c2,
      ZZ(b3),
      ZZ(euler),
      characteristic,
      1,
      0,
      Dict{Symbol, Any}(),
      Dict{Symbol, Any}(),
    )
  end

  function altered_system(S; kwargs...)
    pick(name, value) = haskey(kwargs, name) ? kwargs[name] : value
    return GKMtools._GKMInvariantSystem(
      pick(:H2_rank, S.H2_rank),
      pick(:H2, S.H2),
      pick(:basis_localizations, S.basis_localizations),
      pick(:mu_packed, S.mu_packed),
      pick(:mu_triples, S.mu_triples),
      pick(:cubic_ring, S.cubic_ring),
      pick(:cubic, S.cubic),
      pick(:c1, S.c1),
      pick(:w2, S.w2),
      pick(:p1, S.p1),
      pick(:c2, S.c2),
      pick(:b3, S.b3),
      pick(:euler, S.euler),
      pick(:characteristic_numbers, S.characteristic_numbers),
      pick(:root, S.root),
      pick(:weight_rank, S.weight_rank),
      pick(:diagnostics, Dict{Symbol, Any}()),
      pick(:cache, Dict{Symbol, Any}()),
    )
  end

  rank1 = manual_system([1], [4], [4]; euler=4)
  rank1_c1zero = manual_system([1], [0], [4]; euler=4)
  rank1_p2 = manual_system([1], [4], [28]; euler=4)
  rank1_odd = manual_system([1], [1], [1]; euler=4)
  rank1_zero = manual_system([0], [0], [0])
  rank1_six = manual_system([6], [0], [0])
  rank1_spin_even = manual_system([6], [0], [0])
  rank1_nonspin_even = manual_system([6], [1], [6])
  positive2 = manual_system([1, 0, 0, 1], [2, 2], [0, 0]; euler=6)
  indefinite2 = manual_system([1, 0, 0, -1], [2, 2], [0, 0]; euler=6)
  singular2 = manual_system([0, 0, 0, 0], [0, 0], [0, 0]; euler=0)
  rational_split2 = manual_system([6, 0, 0, 6], [0, 0], [0, 0])
  rational_lines2 = manual_system([0, 6, 6, 0], [0, 0], [0, 0])

  @testset "exports and public graph entry point" begin
    exported = names(GKMtools)
    @test :system_of_invariants_6d in exported
    @test :compare_systems in exported
    @test :pontryagin_class in exported
    @test :first_pontryagin_class in exported
    for internal_name in (
      :_GKMInvariantSystem, :_validate_signed_gkm_graph,
      :_integral_kernel_basis, :_mu_eval, :_verify_system_isomorphism,
      :_transport_invariant_system,
    )
      @test internal_name ∉ exported
    end

    P3 = projective_space(GKM_graph, 3)
    p1 = first_pontryagin_class(P3)
    @test p1 == pontryagin_class(P3, 1)
    @test pontryagin_class(P3, 0) == one(P3.equivariantCohomology)
    @test iszero(pontryagin_class(P3, 4))
    @test_throws ArgumentError pontryagin_class(P3, -1)
    @test p1 == first_chern_class(P3)^2 - 2 * chern_class(P3, 2)
    @test integrate(p1 * first_chern_class(P3), P3) == 16
    S = system_of_invariants_6d(P3)
    @test S.H2_rank == 1
    @test abs(GKMtools._mu_entry(S, 1, 1, 1)) == 1
    @test S.characteristic_numbers == (c1_cubed=ZZ(64), c1_c2=ZZ(24), c3=ZZ(4), p1_c1=ZZ(16))
    @test GKMtools._root_index(P3, 1) == 1
    @test GKMtools._root_index(P3, P3.labels[1]) == 1
    @test_throws ArgumentError GKMtools._root_index(P3, "missing")
    relation_matrix = GKMtools._build_degree_two_relation_matrix(P3, 1)
    @test nrows(relation_matrix) == rank_torus(P3) * length(collect(edges(P3.g)))
    @test ncols(GKMtools._normalized_h2_basis(P3, 1)) == 1
    edge = first(edges(P3.g))
    @test GKMtools._weight_column(P3, edge) ==
      GKMtools._weight_column(P3, src(edge), P3.edge_to_flag_index[edge])

    P2 = projective_space(GKM_graph, 2)
    validation = GKMtools._validate_signed_gkm_graph(P2; check_characteristic_classes=false)
    @test !validation.valid
    @test any(e -> occursin("3-valent", e), validation.errors)
    @test_throws ArgumentError system_of_invariants_6d(P2)
  end

  @testset "manual tensor and localization helpers" begin
    @test GKMtools._symmetric_triples(2) == [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
    @test GKMtools._mu_entry(positive2, 2, 1, 1) == 0
    @test GKMtools._mu_entry(positive2, 2, 2, 2) == 1
    @test_throws ArgumentError GKMtools._mu_entry(positive2, 3, 1, 1)
    @test GKMtools._mu_eval(positive2, [1, 2], [3, 4], [5, 6]) == 63
    @test GKMtools._mu_eval(positive2, matrix(ZZ, 2, 1, [1, 2]), [3, 4], [5, 6]) == 63
    @test_throws ArgumentError GKMtools._mu_eval(positive2, [1], [1], [1])
    @test GKMtools._vector_entry([4, 5], 2) == 5
    @test GKMtools._vector_entry(matrix(ZZ, 2, 1, [4, 5]), 2) == 5

    flattening = GKMtools._mu_flattening(positive2)
    @test flattening == matrix(ZZ, 2, 3, [1, 0, 0, 0, 0, 1])
    @test base_ring(GKMtools._mu_flattening(positive2; base_ring=GF(2))) == GF(2)
    @test GKMtools._contraction_matrix(positive2, [2, 3]) == matrix(ZZ, 2, 2, [2, 0, 0, 3])
    @test GKMtools._multiplicity(1, 1, 1) == 1
    @test GKMtools._multiplicity(1, 1, 2) == 3
    @test GKMtools._multiplicity(1, 2, 3) == 6

    _, cubic = GKMtools._diagonal_cubic(positive2)
    @test evaluate(cubic, [ZZ(2), ZZ(3)]) == 35
    _, determinant = GKMtools._determinant_polynomial(positive2)
    @test evaluate(determinant, [QQ(2), QQ(3)]) == 6

    R, (x,) = polynomial_ring(QQ, ["x"])
    @test GKMtools._extract_integral_constant(R(6) // R(3)) == 2
    @test_throws ArgumentError GKMtools._extract_integral_constant(x // one(R))
    @test_throws ArgumentError GKMtools._extract_integral_constant(R(1) // R(2))
    @test GKMtools._linear_coefficients(3 * x, 1) == matrix(ZZ, 1, 1, [3])
    @test_throws ArgumentError GKMtools._linear_coefficients(x^2, 1)
    @test_throws ArgumentError GKMtools._linear_coefficients((QQ(1) // QQ(2)) * x, 1)

    M = matrix(ZZ, 1, 2, [2, 0])
    K = GKMtools._integral_kernel_basis(M)
    @test M * K == zero_matrix(ZZ, 1, 1)
    @test abs(K[2, 1]) == 1
    @test GKMtools._integral_multiple(matrix(ZZ, 2, 1, [4, 6]), matrix(ZZ, 2, 1, [2, 3]))
    @test !GKMtools._integral_multiple(matrix(ZZ, 2, 1, [1, 1]), matrix(ZZ, 2, 1, [2, 2]))
    @test !GKMtools._integral_multiple(matrix(ZZ, 2, 1, [1, 0]), zero_matrix(ZZ, 2, 1))
  end

  @testset "manual system validation paths" begin
    @test !GKMtools._validate_invariant_system(positive2)
    @test GKMtools._validate_invariant_system(rank1)
    @test GKMtools._validate_invariant_system(rational_split2)
    @test GKMtools._validate_invariant_system(rational_lines2)
    @test GKMtools._realizability_diagnostics(rank1).realizable
    @test !GKMtools._realizability_diagnostics(rank1_odd).realizable
    @test GKMtools._characteristic_numbers(positive2) === positive2.characteristic_numbers
    @test !GKMtools._validate_invariant_system(altered_system(positive2; c1=zero_matrix(ZZ, 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(positive2; w2=zero_matrix(GF(2), 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(positive2; p1=zero_matrix(ZZ, 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(positive2; c2=zero_matrix(ZZ, 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(positive2; mu_packed=zero_matrix(ZZ, 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(positive2; mu_triples=[(1, 1, 1)]))
    @test !GKMtools._validate_invariant_system(altered_system(rank1_odd; w2=zero_matrix(GF(2), 1, 1)))
    @test !GKMtools._validate_invariant_system(altered_system(rank1; c2=rank1.c2 + matrix(ZZ, 1, 1, [1])))

    bad_cubic_number = merge(rank1.characteristic_numbers, (c1_cubed=ZZ(9),))
    bad_c2_number = merge(rank1.characteristic_numbers, (c1_c2=ZZ(5),))
    bad_euler = merge(rank1.characteristic_numbers, (c3=ZZ(5),))
    @test !GKMtools._validate_invariant_system(altered_system(rank1; characteristic_numbers=bad_cubic_number))
    @test !GKMtools._validate_invariant_system(altered_system(rank1; characteristic_numbers=bad_c2_number))
    @test !GKMtools._validate_invariant_system(altered_system(rank1; characteristic_numbers=bad_euler))
  end

  @testset "transport and exact verifier paths" begin
    A = matrix(ZZ, 2, 2, [1, 1, 0, 1])
    transported = GKMtools._transport_invariant_system(rational_split2, A)
    @test GKMtools._verify_system_isomorphism(rational_split2, transported, A)
    @test GKMtools._verify_system_isomorphism(rational_split2, transported, A; preserve_almost_complex=true)
    @test_throws ArgumentError GKMtools._transport_invariant_system(rational_split2, 2 * identity_matrix(ZZ, 2))

    @test !GKMtools._verify_system_isomorphism(rank1, positive2, identity_matrix(ZZ, 1))
    @test !GKMtools._verify_system_isomorphism(positive2, positive2, zero_matrix(ZZ, 1, 1))
    @test !GKMtools._verify_system_isomorphism(positive2, positive2, 2 * identity_matrix(ZZ, 2))
    @test !GKMtools._verify_system_isomorphism(rank1, rank1, matrix(ZZ, 1, 1, [-1]))
    @test !GKMtools._verify_system_isomorphism(rank1, rank1_odd, identity_matrix(ZZ, 1); preserve_almost_complex=false)
    @test !GKMtools._verify_system_isomorphism(rank1, rank1_p2, identity_matrix(ZZ, 1); preserve_almost_complex=false)
    @test !GKMtools._verify_system_isomorphism(rank1, manual_system([3], [2], [0]; euler=4), identity_matrix(ZZ, 1); preserve_almost_complex=false)
    @test !GKMtools._verify_system_isomorphism(
      rank1, altered_system(rank1; b3=ZZ(2)), identity_matrix(ZZ, 1),
    )
    @test !GKMtools._verify_system_isomorphism(
      rank1, altered_system(rank1; euler=ZZ(5)), identity_matrix(ZZ, 1),
    )

  end

  @testset "exact obstruction helpers" begin
    @test GKMtools._content(matrix(ZZ, 1, 3, [0, -6, 9])) == 3
    @test GKMtools._vector_divisibility(matrix(ZZ, 2, 1, [6, -9])) == 3
    @test GKMtools._smith_diagonal(matrix(ZZ, 2, 2, [2, 0, 0, 6])) == ZZ.([2, 6])
    @test GKMtools._tensor_signature(positive2) == GKMtools._tensor_signature(positive2)
    @test GKMtools._contraction_signature(positive2, positive2.c1) == GKMtools._contraction_signature(positive2, positive2.c1)
    @test GKMtools._rational_signature(matrix(ZZ, 3, 3, [1, 0, 0, 0, -1, 0, 0, 0, 0])) == (1, 1, 1)
    @test GKMtools._rational_signature(matrix(ZZ, 2, 2, [0, 2, 2, 0])) == (1, 0, 1)
    @test GKMtools._restricted_tensor_signature(positive2, identity_matrix(ZZ, 2)).content == 1
    @test GKMtools._canonical_sublattice_signatures(positive2; preserve_almost_complex=true) ==
      GKMtools._canonical_sublattice_signatures(positive2; preserve_almost_complex=true)

    @test occursin("H2 ranks", GKMtools._cheap_obstruction(rank1, positive2; preserve_almost_complex=false))
    @test occursin("b3", GKMtools._cheap_obstruction(rank1, altered_system(rank1; b3=ZZ(2)); preserve_almost_complex=false))
    @test occursin("Euler", GKMtools._cheap_obstruction(rank1, altered_system(rank1; euler=ZZ(5)); preserve_almost_complex=false))
    @test occursin("p1 contents", GKMtools._cheap_obstruction(rank1, rank1_p2; preserve_almost_complex=false))
    @test occursin("w2", GKMtools._cheap_obstruction(
      rank1_spin_even, manual_system([6], [1], [0]); preserve_almost_complex=false,
    ))
    @test occursin("tensor", GKMtools._cheap_obstruction(rank1_zero, rank1_six; preserve_almost_complex=false))
    @test occursin("c1 divisibilities", GKMtools._cheap_obstruction(rank1, manual_system([1], [8], [4]; euler=4); preserve_almost_complex=true))
    @test !isnothing(GKMtools._rational_polynomial_obstruction(rational_split2, rational_lines2))
  end

  @testset "finite and bounded search paths" begin
    histogram = GKMtools._finite_field_histogram(positive2, 2; preserve_almost_complex=true)
    @test sum(values(histogram)) == 4
    @test keytype(histogram) == NTuple{8, Int}
    @test all(key -> all(entry -> !(entry isa String), key), keys(histogram))
    @test isnothing(GKMtools._finite_field_signature_obstruction(
      positive2, positive2, 2; preserve_almost_complex=true, point_cap=10,
    ))
    @test !isnothing(GKMtools._finite_field_signature_obstruction(
      rank1_zero, rank1, 2; preserve_almost_complex=false, point_cap=10,
    ))
    @test occursin("contraction", GKMtools._finite_field_signature_obstruction(
      manual_system([1], [0], [0]), rank1_odd, 2;
      preserve_almost_complex=true, point_cap=0,
    ))
    @test sum(values(GKMtools._finite_field_histogram(
      positive2, 3; preserve_almost_complex=false,
    ))) == 9

    status, witnesses, _ = GKMtools._finite_field_isomorphism_search(
      rank1, rank1, 2; preserve_almost_complex=true, node_cap=100,
    )
    @test status == :found
    @test !isempty(witnesses)
    status, _, _ = GKMtools._finite_field_isomorphism_search(
      rank1_zero, rank1, 2; preserve_almost_complex=false, node_cap=100,
    )
    @test status == :none
    status, _, _ = GKMtools._finite_field_isomorphism_search(
      rank1, rank1, 2; preserve_almost_complex=true, node_cap=0,
    )
    @test status == :capped

    F2 = GF(2)
    F3 = GF(3)
    modular_data = [
      (prime=2, solutions=[matrix(F2, 1, 1, [1])], complete=true),
      (prime=3, solutions=[matrix(F3, 1, 1, [2])], complete=true),
    ]
    residues, modulus, complete = GKMtools._crt_residue_classes(modular_data)
    @test complete
    @test modulus == 6
    @test only(residues)[1, 1] == 5

    found = GKMtools._bounded_integral_witness_search(
      rank1, rank1, [1]; preserve_almost_complex=true, node_cap=100,
    )
    @test found.status == :found
    @test found.witness == identity_matrix(ZZ, 1)
    empty_search = GKMtools._bounded_integral_witness_search(
      rank1, rank1, Int[]; preserve_almost_complex=true, node_cap=100,
    )
    @test empty_search.status == :exhausted
    @test isnothing(empty_search.witness)
    capped_search = GKMtools._bounded_integral_witness_search(
      rank1, rank1, [1]; preserve_almost_complex=true, node_cap=0,
    )
    @test capped_search.status == :capped
    @test isnothing(capped_search.witness)
    @test GKMtools._column_matrix([[ZZ(1), ZZ(0)], [ZZ(0), ZZ(1)]]) == identity_matrix(ZZ, 2)
    @test size(GKMtools._column_matrix(Vector{Vector{ZZRingElem}}())) == (0, 0)
    @test GKMtools._partial_tensor_matches(rank1, rank1, [[ZZ(1)]])
    @test !GKMtools._partial_tensor_matches(rank1, rank1, [[ZZ(-1)]])
  end

  @testset "definite-lattice and comparison outcomes" begin
    @test GKMtools._qq_matrix_to_integral(matrix(QQ, 1, 1, [2])) == matrix(ZZ, 1, 1, [2])
    @test isnothing(GKMtools._qq_matrix_to_integral(matrix(QQ, 1, 1, [QQ(1) // QQ(2)])))

    @test GKMtools._definite_contraction_search(
      rank1, rank1; preserve_almost_complex=false,
    )[1] == :skipped
    @test GKMtools._definite_contraction_search(
      singular2, singular2; preserve_almost_complex=true,
    )[1] == :skipped
    @test GKMtools._definite_contraction_search(
      indefinite2, indefinite2; preserve_almost_complex=true,
    )[1] == :skipped
    @test GKMtools._definite_contraction_search(
      rank1, manual_system([-1], [2], [0]; euler=4);
      preserve_almost_complex=true,
    )[1] == :none
    @test GKMtools._definite_contraction_search(
      rank1, manual_system([1], [4], [0]; euler=4);
      preserve_almost_complex=true,
    )[1] == :none
    @test GKMtools._definite_contraction_search(
      rank1, rank1; preserve_almost_complex=true,
    )[1] == :found
    @test GKMtools._definite_contraction_search(
      rank1, rank1_nonspin_even; preserve_almost_complex=true,
    )[1] == :none
    @test GKMtools._definite_contraction_search(
      rank1, rank1; preserve_almost_complex=true,
      group_order_cap=1,
    )[1] == :skipped

    invalid = altered_system(rank1; c1=zero_matrix(ZZ, 2, 1))
    @test_throws ArgumentError compare_systems(invalid, rank1)
    @test_throws ArgumentError compare_systems(rank1, invalid)
    @test_throws ArgumentError compare_systems(rank1_odd, rank1_odd)
    @test compare_systems(rank1_zero, rank1_six; primes=[]).status == :not_equivalent
    @test compare_systems(rank1, rank1_c1zero; primes=[]).status == :equivalent
    @test compare_systems(
      rank1, rank1_c1zero; preserve_almost_complex=true, primes=[],
    ).status == :not_equivalent

    equivalent_definite = compare_systems(rank1, rank1; primes=[])
    @test equivalent_definite.status == :equivalent
    @test GKMtools._verify_system_isomorphism(rank1, rank1, equivalent_definite.witness)
    @test first(equivalent_definite.diagnostics) == (:comparison_mode, :oriented_smooth)
    @test all(key -> !startswith(String(key), "_"), keys(rank1.diagnostics))

    equivalent_bounded = compare_systems(
      singular2, singular2; primes=[], use_definite_contraction=false,
      integral_search_bounds=[1],
    )
    @test equivalent_bounded.status == :equivalent
    equivalent_modular = compare_systems(
      rank1, rank1; primes=[2], use_definite_contraction=false,
      integral_search_bounds=[1],
    )
    @test equivalent_modular.status == :equivalent
    @test any(d -> d[1] == :finite_field_search && d[3] == :found, equivalent_modular.diagnostics)
    capped_modular = compare_systems(
      rank1, rank1; primes=[2], finite_field_point_cap=0,
      finite_field_isomorphism_cap=0, use_definite_contraction=false,
      integral_search_bounds=[1],
    )
    @test capped_modular.status == :equivalent
    @test any(d -> d[1] == :finite_field_search && d[3] == :capped, capped_modular.diagnostics)
    capped_without_residue = compare_systems(
      singular2, singular2; primes=[2], finite_field_point_cap=0,
      finite_field_isomorphism_cap=0, use_definite_contraction=false,
      integral_search_bounds=Int[],
    )
    @test capped_without_residue.status == :unknown
    integral_diagnostic = only(filter(
      d -> d[1] == :integral_search,
      capped_without_residue.diagnostics,
    ))
    @test length(integral_diagnostic[5]) == 1
    @test only(integral_diagnostic[5])[1] == :unconstrained
    @test any(d -> d[1] == :integral_search, compare_systems(
      singular2, singular2; primes=[], use_definite_contraction=false,
      integral_search_bounds=Int[],
    ).diagnostics)

    unknown = compare_systems(
      singular2, singular2; primes=[], use_definite_contraction=false,
      integral_search_bounds=Int[],
    )
    @test unknown.status == :unknown
    @test isnothing(unknown.witness)
    @test any(d -> d[1] == :integral_search && d[2] == :exhausted_bounds, unknown.diagnostics)

    @test_throws ArgumentError GKMtools._SystemComparisonResult(:equivalent, nothing, nothing, Any[])
    @test_throws ArgumentError GKMtools._SystemComparisonResult(:not_equivalent, nothing, nothing, Any[])

    io = IOBuffer()
    show(io, unknown)
    @test occursin("unknown", String(take!(io)))
  end

  @testset "known six-manifold systems and transported bases" begin
    product_graph(G1, G2) = *(
      G1, G2;
      calculateCurveClasses=false,
      calculateConnection=false,
    )
    P1 = projective_space(GKM_graph, 1)
    fixtures = (
      (
        projective_space(GKM_graph, 3),
        manual_system([1], [4], [4]; euler=4),
      ),
      (
        product_graph(product_graph(P1, P1), P1),
        manual_system([0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [2, 2, 2], [0, 0, 0]; euler=8),
      ),
      (
        product_graph(P1, projective_space(GKM_graph, 2)),
        manual_system([0, 0, 1, 0], [2, 3], [3, 0]; euler=6),
      ),
    )

    computed_systems = GKMtools._GKMInvariantSystem[]
    for (graph, expected) in fixtures
      computed = system_of_invariants_6d(graph)
      push!(computed_systems, computed)
      witness = -identity_matrix(ZZ, expected.H2_rank)
      @test GKMtools._verify_system_isomorphism(
        expected, computed, witness; preserve_almost_complex=true,
      )
      @test computed.b3 == 0
      @test computed.diagnostics[:realizability].realizable
    end

    root_graph = fixtures[3][1]
    root_system = computed_systems[3]
    for root in vertices(root_graph.g)
      rerooted = system_of_invariants_6d(root_graph; root)
      comparison = compare_systems(
        root_system, rerooted;
        preserve_almost_complex=true,
        primes=[2],
        integral_search_bounds=[1, 2],
      )
      @test comparison.status == :equivalent
      @test GKMtools._verify_system_isomorphism(
        root_system, rerooted, comparison.witness; preserve_almost_complex=true,
      )
    end

    base = computed_systems[2]
    transports = (
      matrix(ZZ, 3, 3, [1, 1, 0, 0, 1, 0, 0, 0, 1]),
      matrix(ZZ, 3, 3, [1, 5, 0, 0, 1, 0, 0, 0, 1]),
      matrix(ZZ, 3, 3, [1, 2, 3, 0, 1, 4, 0, 0, 1]),
    )
    for A in transports
      transported = GKMtools._transport_invariant_system(base, A)
      @test GKMtools._verify_system_isomorphism(
        base, transported, A; preserve_almost_complex=true,
      )
      result = compare_systems(
        base, transported;
        preserve_almost_complex=true,
        primes=[2],
        integral_search_bounds=[1, 2],
      )
      @test result.status != :not_equivalent
      if result.status == :equivalent
        @test GKMtools._verify_system_isomorphism(
          base, transported, result.witness; preserve_almost_complex=true,
        )
      end
    end

    function rebuild_graph(G; permutation=collect(1:n_vertices(G.g)), weight_change=nothing)
      n = n_vertices(G.g)
      k = rank_torus(G)
      isnothing(weight_change) && (weight_change = identity_matrix(ZZ, k))
      inverse_permutation = invperm(permutation)
      graph = Graph{Undirected}(n)
      for edge in edges(G.g)
        add_edge!(graph, permutation[src(edge)], permutation[dst(edge)])
      end
      lattice = free_module(ZZ, k)
      lattice_basis = gens(lattice)
      weights = Dict{Edge, elem_type(lattice)}()
      for edge in edges(graph)
        old_source = inverse_permutation[src(edge)]
        old_target = inverse_permutation[dst(edge)]
        old_weight = GKMtools._weight_column(G, Edge(old_source, old_target))
        new_weight = weight_change * old_weight
        weights[edge] = sum(
          new_weight[i, 1] * lattice_basis[i] for i in 1:k; init=zero(lattice),
        )
      end
      labels = ["v$i" for i in 1:n]
      return gkm_graph(graph, labels, lattice, weights)
    end

    original_graph = fixtures[1][1]
    original_system = computed_systems[1]
    reversed_vertices = rebuild_graph(original_graph; permutation=reverse(collect(1:4)))
    k = rank_torus(original_graph)
    shear = identity_matrix(ZZ, k)
    k >= 2 && (shear[1, 2] = 3)
    changed_weights = rebuild_graph(original_graph; weight_change=shear)
    for rebuilt in (reversed_vertices, changed_weights)
      rebuilt_system = system_of_invariants_6d(rebuilt)
      result = compare_systems(
        original_system, rebuilt_system;
        preserve_almost_complex=true,
        primes=[2],
        integral_search_bounds=[1, 2],
      )
      @test result.status == :equivalent
      @test GKMtools._verify_system_isomorphism(
        original_system, rebuilt_system, result.witness; preserve_almost_complex=true,
      )
    end

    unchecked = system_of_invariants_6d(original_graph; check=false)
    @test occursin("skipped", only(unchecked.diagnostics[:validation].warnings))

    zero_weight = deepcopy(original_graph)
    zero_weight.weights_at_vertex[1][1] = zero(zero_weight.M)
    zero_validation = GKMtools._validate_signed_gkm_graph(
      zero_weight; check_characteristic_classes=false,
    )
    @test !zero_validation.valid

    endpoint_mismatch = deepcopy(original_graph)
    edge = first(edges(endpoint_mismatch.g))
    reverse_index = endpoint_mismatch.edge_to_flag_index[reverse(edge)]
    endpoint_mismatch.weights_at_vertex[dst(edge)][reverse_index] += gens(endpoint_mismatch.M)[1]
    endpoint_validation = GKMtools._validate_signed_gkm_graph(
      endpoint_mismatch; check_characteristic_classes=false,
    )
    @test !endpoint_validation.valid

    dependent = deepcopy(original_graph)
    dependent.weights_at_vertex[1][2] = dependent.weights_at_vertex[1][1]
    dependent_validation = GKMtools._validate_signed_gkm_graph(
      dependent; check_characteristic_classes=false,
    )
    @test !dependent_validation.valid

    characteristic_failure_found = false
    for changed_edge in edges(original_graph.g), basis_index in 1:k
      trial = rebuild_graph(original_graph)
      source_index = trial.edge_to_flag_index[changed_edge]
      target_index = trial.edge_to_flag_index[reverse(changed_edge)]
      changed_weight = trial.weights_at_vertex[src(changed_edge)][source_index] +
        gens(trial.M)[basis_index]
      trial.weights_at_vertex[src(changed_edge)][source_index] = changed_weight
      trial.weights_at_vertex[dst(changed_edge)][target_index] = -changed_weight
      isvalid(trial; printDiagnostics=false) || continue
      trial_validation = GKMtools._validate_signed_gkm_graph(trial)
      if !trial_validation.valid && any(
        error -> occursin("c1", error) || occursin("p1", error),
        trial_validation.errors,
      )
        characteristic_failure_found = true
        break
      end
    end
    @test characteristic_failure_found

    scaled_weights = rebuild_graph(original_graph; weight_change=2 * identity_matrix(ZZ, k))
    scaled_validation = GKMtools._validate_signed_gkm_graph(scaled_weights)
    @test scaled_validation.valid
    @test any(warning -> occursin("non-primitive", warning), scaled_validation.warnings)

    rational_result = compare_systems(rational_split2, rational_lines2; primes=[])
    @test rational_result.status == :not_equivalent
    @test rational_result.obstruction[1] == :rational_factorization

    io = IOBuffer()
    show(io, computed_systems[1])
    displayed = String(take!(io))
    @test occursin("rank(H^2) = 1", displayed)
    @test !occursin("basis_localizations", displayed)
  end
end
