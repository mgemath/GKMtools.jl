@testset "Batched small quantum product" begin
    G = projective_space(GKMGraph, 1)
    h = first_chern_class(G) / 2
    basis = [h^0, h]
    line = curve_class(G, first(GKMtools.edges(G)))
    classes = [h^0, h, 2h, h]

    for equivariant in (false, true)
        for beta in (zero(parent(line)), line, nothing)
            matrices = matrix_small_quantum_product(G, classes;
                basis, beta, equivariant, show_bar=false)
            @test length(matrices) == length(classes)
            @test matrices == [matrix_small_quantum_product(G, c;
                basis, beta, equivariant, show_bar=false) for c in classes]
            @test matrices[3] == 2matrices[2]
            @test matrices[2] == matrices[4]
            @test matrices[2] !== matrices[4]
        end
    end

    classical = matrix_small_quantum_product(G, h; basis, beta=zero(parent(line)), show_bar=false)
    quantum = matrix_small_quantum_product(G, h; basis, beta=line, show_bar=false)
    @test classical == matrix(QQ, [0 0; 1 0])
    @test quantum == matrix(QQ, [0 1; 0 0])
    @test matrix_small_quantum_product(G, h; basis, show_bar=false) == classical + quantum
    @test isempty(matrix_small_quantum_product(G, empty(classes); show_bar=false))
    other = projective_space(GKMGraph, 1)
    @test_throws ArgumentError matrix_small_quantum_product(G, [first_chern_class(other)]; show_bar=false)
end

@testset "Unmarked matrix invariants match marked integration" begin
    G = projective_space(GKMGraph, 2)
    h = first_chern_class(G) / 3
    basis = [h^0, h, h^2]
    classes = [h, h^2]
    line = curve_class(G, first(GKMtools.edges(G)))
    indices = GKMtools.get_symmetric_indices(length(basis))
    marked = [ev(1,c)*ev(2,basis[j])*ev(3,basis[k])
              for c in classes for (j,k) in indices]
    for equivariant in (false, true)
        @test GKMtools._invariants(classes, basis, indices, G, line, false, equivariant) ==
            gromov_witten(G, line, 3, marked; show_bar=false, fast_mode=!equivariant)
    end
end
