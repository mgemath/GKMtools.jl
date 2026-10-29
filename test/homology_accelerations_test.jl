# Run independently with: julia --project -e 'using Test, Oscar, GKMtools; include("test/homology_accelerations_test.jl")'
function full_homology_reference(G)
    data = core(G)
    es, indices = GKMtools._edge_list_and_indices(data)
    E = free_module(ZZ, length(es))
    R = GKMtools._homology_relations(data, es, indices; project=false)
    relations, _ = sub(E, [E([R[i,j] for j in eachindex(es)]) for i in 1:nrows(R)])
    Q, q = quo(E, relations)
    H, p = GKMtools._remove_homology_torsion(Q)
    return H, compose(q,p)
end

@testset "Homology accelerations" begin
    examples = [projective_space(GKMGraph, 1), projective_space(GKMGraph, 3),
        projective_space(GKMGraph,1)*projective_space(GKMGraph,2),
        flag_variety(GKMGraph, [1,1,1]),
        gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety,5)),
        total_space(tangent_bundle(projective_space(GKMGraph,1))),
        total_space(dual(tangent_bundle(projective_space(GKMGraph,2)))),
        total_space(dual(tangent_bundle(projective_space(GKMGraph,1)*projective_space(GKMGraph,1))))]
    for G in examples
        H, q = full_homology_reference(G)
        actual = GKM_second_homology(G)
        @test rank(actual.H2) == rank(H)
        # Equality of the saturated relation lattices identifies every edge
        # class, integrally and independently of the chosen homology bases.
        @test hnf(matrix(last(kernel(q)))) == hnf(matrix(last(kernel(actual.quotient))))
        @test all(chern_number(G,curve_class(G,e)) == GKMtools._homology_chern_number(core(G),e)
                  for e in GKMtools.edges(G))
    end
    @test rank(examples[end-1].H2.H2) == 1
    @test rank(examples[end].H2.H2) == 2
    @test all(iszero(chern_number(examples[end-1],curve_class(examples[end-1],e)))
              for e in GKMtools.edges(examples[end-1]))
    @test !GKMtools._homology_use_compact_formula(core(examples[end-1]))
    @test_throws ArgumentError betti_numbers(examples[end-1])

    G = examples[3]
    P = GKMtools._homology_rank_two_projection(core(G))
    @test P == GKMtools._homology_rank_two_projection(core(G))
    @test det(P[:,1:2]) == 1
    for v in GKMtools.vertices(G)
        ws = [P*matrix(QQ,rank_torus(G),1,[weight(G,v,i)[j] for j in 1:rank_torus(G)]) for i in eachindex(flags(G,v))]
        @test all(rank(hcat(ws[i],ws[j])) == 2 for i in eachindex(ws) for j in (i+1):length(ws))
    end
    @test GKMtools._homology_rank_two_projection(core(examples[1])) == identity_matrix(QQ,2)
end

# Synthetic rank-one data preserves the supplied coordinates, including large
# weights, so it tests the enumeration boundary independently of construction.
function partition_test_data(weights)
    E = free_module(ZZ, length(weights))
    H = free_module(ZZ, 1)
    q = ModuleHomomorphism(E, H, [ZZ(w)*gen(H,1) for w in weights])
    C = cone_from_inequalities(matrix(QQ,1,1,[-1]))
    z = free_module(ZZ,1)
    c = ModuleHomomorphism(H,z,[zero(z)])
    es = [Edge(i,i+1) for i in eachindex(weights)]
    data = GKM_H2(E,H,Dict(es .=> eachindex(weights)),q,C,QQFieldElem[1//minimum(abs.(ZZ.(weights)))],c)
    return data, es, gen(H,1)
end

@testset "Bounded rank-one partition arithmetic" begin
    # Compare the actual production routine with exhaustive exact arithmetic.
    for k in 1:3, ws in Iterators.product(ntuple(_ -> 1:3, k)...)
        weights = collect(ws)
        data, es, b = partition_test_data(weights)
        for d in 1:10
            expected = Set(collect(m) for m in Iterators.product(ntuple(_ -> 1:d, k)...)
                if sum(ZZ(weights[i])*m[i] for i in 1:k) == d)
            @test GKMtools._multiplicities(data,es,d*b) == expected
        end
    end
    M = typemax(Int)
    for (weights, d, expected) in [([M,M,3], 1, Set{Vector{Int}}()),
                                  ([M-1,1], M, Set([[1,1]])),
                                  ([1], M, Set([[M]])),
                                  ([-2,-3], -7, Set([[2,1]]))]
        data, es, b = partition_test_data(weights)
        @test GKMtools._multiplicities(data,es,ZZ(d)*b) == expected
    end

    n = ZZ(2)^100
    data, es, b = partition_test_data([n,n+1])
    @test_throws OverflowError GKMtools._multiplicities(data,es,(2n+1)*b)
    @test_throws r"curve degree.*Int range" GKMtools._multiplicities(data,es,(2n+1)*b)
    @test_throws r"edge weight.*Int range" GKMtools._multiplicities(data,es,b)
    @test_throws r"edge weight.*Int range" GKMtools._gw_max_n_edges(data,b)
    @test GKMtools._max_n_edges(data,n^2*b) == n # Homology bound stays exact.

    data, es, b = partition_test_data([1])
    @test_throws OverflowError GKMtools._multiplicities(data,es,(ZZ(M)+1)*b)
    @test_throws OverflowError GKMtools._multiplicities(data,es,ZZ(typemin(Int))*b)
    @test GKMtools._gw_max_n_edges(data,(ZZ(M)-1)*b) == M-1
    @test_throws r"graph-size bound.*Int range" GKMtools._gw_max_n_edges(data,ZZ(M)*b)
    @test_throws OverflowError GKMtools._gw_max_n_edges(data,(ZZ(M)-1)*b; extra=2)

    # Exercise both GW entry paths and the fast-mode degree precheck. Large
    # edge weights must fail even when the normalized edge bound is tiny.
    P1 = projective_space(GKMGraph,1)
    data, es, b = partition_test_data([n])
    large_weight_graph = typeof(P1)(core(P1),GKMtools.connection(P1),get_cohomology(P1),data,nothing)
    for G in (P1, large_weight_graph)
        beta = G === P1 ? n*curve_class(P1,first(GKMtools.edges(P1))) : b
        for fast in (false,true), genus in (0,1)
            @test_throws OverflowError gromov_witten(G,beta,0,class_one();show_bar=false,fast_mode=fast,g=genus)
        end
        @test_throws OverflowError GKMtools._must_return_zero(G,beta,0,[class_one()])
    end
end

@testset "GW enumeration with bounded Int partitions" begin
    L = total_space(vector_bundle_O(1,[-1,-1];big_torus=true))
    b = curve_class(L,first(GKMtools.edges(L)))
    @test gromov_witten(L,3b,0,class_one();show_bar=false,fast_mode=true) == QQ(1,27)
    F = gkm_3d_twisted_flag()
    b = curve_class(F,"3","4")
    @test gromov_witten(F,b,0,class_one();show_bar=false,g=1) == QQ(1,12)
end
