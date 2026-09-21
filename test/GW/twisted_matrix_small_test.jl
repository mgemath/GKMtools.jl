@testset "Batched twisted small quantum product" begin
    V = vector_bundle_O(2, [1]; enlarge_torus=false)
    G = baseof(V)
    h = first_chern_class(V)
    basis = [h^0, h]
    line = curve_class(G, first(GKMtools.edges(G)))
    classes = [h^0, h, 2h, h]
    indices = GKMtools.get_symmetric_indices(length(basis))

    for c in [classes; [h^2, h^0 + h]]
        insertions = [ev(1,c)*ev(2,basis[j])*ev(3,basis[k])*virtual_zero_section(V)
                      for (j,k) in indices]
        shortcut = GKMtools._twisted_c1_known_zero(V, line, c, basis, indices)
        @test isnothing(shortcut) || shortcut == GKMtools._must_return_zero(G, line, 3, insertions)
    end

    for equivariant in (false, true)
        for beta in (zero(parent(line)), line)
            matrices = twisted_matrix_small_quantum_product(V, classes;
                basis, beta, equivariant, show_bar=false)
            @test matrices == [twisted_matrix_small_quantum_product(V, c;
                basis, beta, equivariant, show_bar=false) for c in classes]
            @test matrices[3] == 2matrices[2]
            @test matrices[2] == matrices[4]
            @test matrices[2] !== matrices[4]
        end
    end
    classical = twisted_matrix_small_quantum_product(V, classes;
        basis, beta=zero(parent(line)), show_bar=false)
    quantum = twisted_matrix_small_quantum_product(V, classes;
        basis, beta=line, show_bar=false)
    @test classical[2] == matrix(QQ, [0 0; 1 0])
    @test quantum[2] == matrix(QQ, [0 1; 0 0])
    @test twisted_matrix_small_quantum_product(V, classes; basis, show_bar=false) == classical + quantum
    @test all(iszero, twisted_matrix_small_quantum_product(V, classes;
        basis, beta=2line, show_bar=false))
    mixed = [h, h^0 + h]
    @test twisted_matrix_small_quantum_product(V, mixed; basis, beta=line, show_bar=false) ==
        [twisted_matrix_small_quantum_product(V, c; basis, beta=line, show_bar=false) for c in mixed]
    @test isempty(twisted_matrix_small_quantum_product(V, empty(classes); show_bar=false))
    other = projective_space(GKMGraph, 2)
    @test_throws ArgumentError twisted_matrix_small_quantum_product(V,
        [first_chern_class(other)]; show_bar=false)
end

@testset "Unmarked twisted invariants match marked integration" begin
    V = vector_bundle_O(2, [1]; enlarge_torus=false)
    G = baseof(V)
    h = first_chern_class(V)
    basis = [h^0, h]
    classes = [h, h^2, h^0 + h]
    line = curve_class(G, first(GKMtools.edges(G)))
    indices = GKMtools.get_symmetric_indices(length(basis))
    twist = virtual_zero_section(V)
    marked = [ev(1,c)*ev(2,basis[j])*ev(3,basis[k])*twist
              for c in classes for (j,k) in indices]
    for equivariant in (false, true)
        @test GKMtools._twisted_invariants(classes, basis, indices, G, line, false, equivariant, twist, V) ==
            gromov_witten(G, line, 3, marked; show_bar=false, fast_mode=!equivariant)
    end
end
