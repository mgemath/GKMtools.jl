using Test
using Oscar
using GKMtools

@testset "exceptional-divisor connections" begin
  affine_sp = affine_space(NormalToricVariety, 4)
  ambient = gkm_graph_of_toric(affine_sp)
  center = subgraph_from_vertices(ambient, [1])

  ordinary_blowup = blow_up(center)
  @test ordinary_blowup isa GKMSubgraph
  exceptional = subgraph(ordinary_blowup)
  @test !isempty(GKMtools.transport(GKMtools.connection(exceptional)))
  @test GKMtools.is_valid(exceptional, GKMtools.connection(exceptional))
  @test GKMtools.is_valid(
    ambient_graph(ordinary_blowup),
    GKMtools.connection(ambient_graph(ordinary_blowup)),
  )

  weighted_blowup = blow_up(center, [1, 2, 3, 4])
  @test weighted_blowup isa OrbifoldGKMSubgraph
  @test ambient_graph(weighted_blowup) isa GKMtools.OrbifoldGKMGraph
  @test subgraph(weighted_blowup) isa GKMtools.OrbifoldGKMGraph
  weighted_exceptional = subgraph(weighted_blowup)
  @test base_ring(GKMtools.lattice(weighted_exceptional)) == ZZ
  @test !isempty(GKMtools.transport(GKMtools.connection(weighted_exceptional)))
  @test GKMtools.is_valid(
    weighted_exceptional,
    GKMtools.connection(weighted_exceptional),
  )
  weighted_ambient = ambient_graph(weighted_blowup)
  @test !isempty(GKMtools.transport(GKMtools.connection(weighted_ambient)))
  @test GKMtools.is_valid(
    weighted_ambient,
    GKMtools.connection(weighted_ambient),
  )

  exceptional = subgraph(weighted_blowup)
  weighted_projective = gkm_graph_of_orbifold_toric(
    stacky_weighted_projective_space_fan([1, 2, 3, 4]),
  )
  isotropy_order(I) = prod(I.isotropy_group; init=1)
  vertex_permutation = [
    findfirst(
      v -> isotropy_order(weighted_projective.vertex_isotropy[v]) ==
           isotropy_order(exceptional.vertex_isotropy[u]),
      vertices(weighted_projective),
    )
    for u in vertices(exceptional)
  ]
  @test vertex_permutation == [1, 4, 3, 2]
  coordinates(weight) = [weight[i] for i in 1:rank(parent(weight))]
  permute_coordinates(weight) = coordinates(weight)[[2, 3, 4, 1]]
  for u in vertices(exceptional)
    v = vertex_permutation[u]
    @test sort(permute_coordinates.(getfield.(exceptional.core.flags[u], :weight))) ==
          sort(coordinates.(getfield.(weighted_projective.core.flags[v], :weight)))
    @test sort(vec(collect(exceptional.vertex_isotropy[u].tangent_rep))) ==
          sort(vec(collect(weighted_projective.vertex_isotropy[v].tangent_rep)))
    @test sort(isotropy_order.(exceptional.flag_isotropy[u])) ==
          sort(isotropy_order.(weighted_projective.flag_isotropy[v]))
    @test sort([(nrows(I.embedding), ncols(I.embedding), Tuple(vec(collect(I.embedding)))) for I in exceptional.flag_isotropy[u]]) ==
          sort([(nrows(I.embedding), ncols(I.embedding), Tuple(vec(collect(I.embedding)))) for I in weighted_projective.flag_isotropy[v]])
  end

  orbifold_center = subgraph_from_vertices(weighted_projective, [1])
  orbifold_blowup = @test_nowarn blow_up(orbifold_center)
  @test orbifold_blowup isa OrbifoldGKMSubgraph
  orbifold_blowup_ambient = ambient_graph(orbifold_blowup)
  @test GKMtools.is_valid(
    orbifold_blowup_ambient, GKMtools.connection(orbifold_blowup_ambient),
  )
  @test sort([
    isotropy_order(I) for I in orbifold_blowup_ambient.vertex_isotropy
    if !isempty(I.isotropy_group)
  ]) == [2, 3, 4]
end
