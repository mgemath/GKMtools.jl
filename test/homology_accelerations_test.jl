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

@testset "Exact rank-one partition arithmetic" begin
    E = free_module(ZZ,2)
    n = ZZ(2)^100
    H, q = GKMtools._primitive_homology_quotient(E,[n,n+1])
    C = cone_from_inequalities(matrix(QQ,1,1,[-1]))
    z = free_module(ZZ,1)
    c = ModuleHomomorphism(H,z,[zero(z)])
    es = [Edge(1,2),Edge(2,3)]
    data = GKM_H2(E,H,Dict(es .=> [1,2]),q,C,QQFieldElem[1//n],c)
    @test GKMtools._multiplicities(data,es,(2n+1)*gen(H,1)) == Set([[1,1]])
    @test GKMtools._multiplicities(data,es,(2n)*gen(H,1)) == Set{Vector{Int}}()
    @test GKMtools._multiplicities(data,[es[1]],3n*gen(H,1)) == Set([[3]])
    @test_throws InexactError GKMtools._multiplicities(data,[es[1]],n^2*gen(H,1))
    @test GKMtools._max_n_edges(data,n^2*gen(H,1)) == n
end
