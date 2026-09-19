@testset "Twisted c1 dimension shortcut" begin
    V = vector_bundle_O(7, [2, 2, 2]; enlarge_torus=false)
    G = baseof(V)
    H = polynomial_class(G, collect(gens(equivariant_coefficient_ring(G))))
    basis = [H^i for i in 0:4]
    class = first_chern_class(G) - first_chern_class(V)
    indices = [(j,k) for j in 1:5 for k in j:5]
    line = curve_class(G, Edge(1,2))
    for d in 1:7
        @test GKMtools._twisted_c1_known_zero(V, d*line, class, basis, indices) ==
            [i+j-3-2*d != 0 for i in 0:4 for j in i:4]
    end
    @test iszero(twisted_c1_matrix(V, 4*line; basis, fast_mode=true))
    @test isnothing(GKMtools._twisted_c1_known_zero(V, line, class, [one(H)+H], [(1,1)]))
    W = vector_bundle_O(2, [1]; enlarge_torus=false)
    X = baseof(W)
    h = polynomial_class(X, collect(gens(equivariant_coefficient_ring(X))))
    b = [h^i for i in 0:1]
    c = first_chern_class(X)-first_chern_class(W)
    beta = curve_class(X, Edge(1,2))
    inds = [(j,k) for j in 1:2 for k in j:2]
    insertions = [ev(1,c)*ev(2,b[j])*ev(3,b[k])*virtual_zero_section(W) for (j,k) in inds]
    @test GKMtools._twisted_c1_known_zero(W,beta,c,b,inds) == GKMtools._must_return_zero(X,beta,3,insertions)
end
