_product_vertex_index(v1::Int, v2::Int, n1::Int) = v1 + (v2 - 1) * n1

function _product_lattice_maps(G1::AbstractGKMGraph{R}, G2::AbstractGKMGraph{R}) where {R}
  M1, M2 = lattice(G1), lattice(G2)
  M = free_module(base_ring(M1), rank(M1) + rank(M2))
  f1 = hom(M1, M, gens(M)[1:rank(M1)])
  f2 = hom(M2, M, gens(M)[rank(M1) + 1:end])
  return M, f1, f2
end

function _product_graph_and_edges(G1, G2)
  n1, n2 = num_vertices(G1), num_vertices(G2)
  g = Graph{Undirected}(n1 * n2)
  origins = Dict{Edge,Tuple{Int,Edge}}()
  for e in edges(G1), v2 in 1:n2
    u = _product_vertex_index(src(e), v2, n1)
    v = _product_vertex_index(dst(e), v2, n1)
    add_edge!(g, u, v)
    E = Edge(u, v) in edges(g) ? Edge(u, v) : Edge(v, u)
    original = src(E) == u ? e : reverse(e)
    origins[E], origins[reverse(E)] = (1, original), (1, reverse(original))
  end
  for e in edges(G2), v1 in 1:n1
    u = _product_vertex_index(v1, src(e), n1)
    v = _product_vertex_index(v1, dst(e), n1)
    add_edge!(g, u, v)
    E = Edge(u, v) in edges(g) ? Edge(u, v) : Edge(v, u)
    original = src(E) == u ? e : reverse(e)
    origins[E], origins[reverse(E)] = (2, original), (2, reverse(original))
  end
  return g, origins
end

function _product_labels(G1, G2)
  n1, n2 = num_vertices(G1), num_vertices(G2)
  V1, V2 = eltype(labels(G1)), eltype(labels(G2))
  result = Vector{ProductVertex{V1,V2}}(undef, n1 * n2)
  for v2 in 1:n2, v1 in 1:n1
    a, b = labels(G1)[v1], labels(G2)[v2]
    result[_product_vertex_index(v1, v2, n1)] =
      ProductVertex(a, b)
  end
  return result
end

function _product_edge_flags(G1, G2, g, origins)
  val1 = valency(G1)
  result = Dict{Edge,Tuple{Int,Int}}()
  for E in edges(g)
    factor, e = origins[E]
    if factor == 1
      result[E] = _edge_flag_indices(G1, e)
    else
      i, j = _edge_flag_indices(G2, e)
      result[E] = (val1 + i, val1 + j)
    end
  end
  return result
end

function _product_connection(G1::AbstractGKMGraph, G2::AbstractGKMGraph,
                             origins, val1::Int, val2::Int, ::Type{C}) where {C}
  trans = Dict{Edge,Vector{Int}}()
  coeffs = Dict{Edge,Vector{C}}()
  for (E, (factor, e)) in origins
    if factor == 1
      trans[E] = vcat(transport(G1)[e], collect(val1 + 1:val1 + val2))
      coeffs[E] = vcat(C.(coefficients(G1)[e]), zeros(C, val2))
    else
      trans[E] = vcat(collect(1:val1), val1 .+ transport(G2)[e])
      coeffs[E] = vcat(zeros(C, val1), C.(coefficients(G2)[e]))
    end
  end
  return Connection{C}(trans, coeffs, "Product($(connection_type(G1)), $(connection_type(G2)))")
end

function _block_diagonal(A::ZZMatrix, B::ZZMatrix)
  result = zero_matrix(ZZ, nrows(A) + nrows(B), ncols(A) + ncols(B))
  for i in 1:nrows(A), j in 1:ncols(A)
    result[i, j] = A[i, j]
  end
  for i in 1:nrows(B), j in 1:ncols(B)
    result[nrows(A) + i, ncols(A) + j] = B[i, j]
  end
  return result
end

@doc raw"""
    *(G1::GKMGraph{R}, G2::GKMGraph{R}) where {R} -> GKMGraph{R}


Construct the Cartesian product of two smooth GKM graphs, or of two
orbifold GKM graphs. The character lattices are combined by direct sum.

# Examples
```jldoctest
julia> G = generalized_gkm_flag(root_system(:A, 1))
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)
Birkhoff-Grothendieck connection for GKM graph with 2 nodes and valency 1

julia> G*G
GKM graph with 4 nodes, valency 2 and axial function:
s1,id -> id,id => (-1, 1, 0, 0)
id,s1 -> id,id => (0, 0, -1, 1)
s1,s1 -> s1,id => (0, 0, -1, 1)
s1,s1 -> id,s1 => (-1, 1, 0, 0)
Product(Birkhoff-Grothendieck, Birkhoff-Grothendieck) connection for GKM graph with 4 nodes and valency 2

julia> W = stacky_weighted_projective_space_fan([1, 2, 4]);

julia> Wp = gkm_graph_of_orbifold_toric(W);

julia> Wp*Wp
Orbifold GKM graph with 9 nodes, valency 4 and axial function:
2,1 -> 1,1 => (0, -1, 2, 0, 0, 0)
3,1 -> 1,1 => (-1, 0, 1, 0, 0, 0)
3,1 -> 2,1 => (-2, 1, 0, 0, 0, 0)
1,2 -> 1,1 => (0, 0, 0, 0, -1, 2)
2,2 -> 2,1 => (0, 0, 0, 0, -1, 2)
2,2 -> 1,2 => (0, -1, 2, 0, 0, 0)
3,2 -> 3,1 => (0, 0, 0, 0, -1, 2)
3,2 -> 1,2 => (-1, 0, 1, 0, 0, 0)
3,2 -> 2,2 => (-2, 1, 0, 0, 0, 0)
1,3 -> 1,1 => (0, 0, 0, -1, 0, 1)
1,3 -> 1,2 => (0, 0, 0, -2, 1, 0)
2,3 -> 2,1 => (0, 0, 0, -1, 0, 1)
2,3 -> 2,2 => (0, 0, 0, -2, 1, 0)
2,3 -> 1,3 => (0, -1, 2, 0, 0, 0)
3,3 -> 3,1 => (0, 0, 0, -1, 0, 1)
3,3 -> 3,2 => (0, 0, 0, -2, 1, 0)
3,3 -> 1,3 => (-1, 0, 1, 0, 0, 0)
3,3 -> 2,3 => (-2, 1, 0, 0, 0, 0)
Vertex Isotropy:
2,1 => Orbifold Vertex Isotropy, Group Structure (cyclic factors): [2]
  Tangent Representation: 
[1   1   0   0]
1,2 => Orbifold Vertex Isotropy, Group Structure (cyclic factors): [2]
  Tangent Representation: 
[0   0   1   1]
2,2 => Orbifold Vertex Isotropy, Group Structure (cyclic factors): [2, 2]
  Tangent Representation: 
[1   1   0   0]
[0   0   1   1]
3,2 => Orbifold Vertex Isotropy, Group Structure (cyclic factors): [2]
  Tangent Representation: 
[0   0   1   1]
2,3 => Orbifold Vertex Isotropy, Group Structure (cyclic factors): [2]
  Tangent Representation: 
[1   1   0   0]
Flag Isotropy:
2,1.3 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,1.4 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
1,2.1 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
1,2.2 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,2.1 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,2.2 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,2.3 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,2.4 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
3,2.1 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
3,2.2 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,3.3 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
2,3.4 => Orbifold Flag Isotropy, Group Structure (cyclic factors): [2]
  Embedding Matrix: 
[1]
Product(Algorithmic, Algorithmic) connection for GKM graph with 9 nodes and valency 4
```
"""
function *(G1::GKMGraph{R}, G2::GKMGraph{R}) where {R}
  n1, n2 = num_vertices(G1), num_vertices(G2)
  val1, val2 = valency(G1), valency(G2)
  M, f1, f2 = _product_lattice_maps(G1, G2)
  g, origins = _product_graph_and_edges(G1, G2)
  product_flags = Vector{Vector{FlagWeight{R}}}(undef, n1 * n2)
  for v2 in 1:n2, v1 in 1:n1
    V = _product_vertex_index(v1, v2, n1)
    product_flags[V] = vcat(
      [FlagWeight{R}(f1(weight(G1, v1, i))) for i in 1:val1],
      [FlagWeight{R}(f2(weight(G2, v2, i))) for i in 1:val2],
    )
  end
  product_labels = _product_labels(G1, G2)
  core = GKMCombinatorialData{R,eltype(product_labels),FlagWeight{R}}(
    g, M, product_labels, product_flags, _product_edge_flags(G1, G2, g, origins),
  )
  con = _product_connection(G1, G2, origins, val1, val2, R)
  return GKMGraph{R,eltype(product_labels),FlagWeight{R}}(
    core, con, create_cohomology(rank(M), n1 * n2), nothing,
  )
end

function Base.:*(G1::OrbifoldGKMGraph{R}, G2::OrbifoldGKMGraph{R}) where {R}
  n1, n2 = num_vertices(G1), num_vertices(G2)
  val1, val2 = valency(G1), valency(G2)
  M, f1, f2 = _product_lattice_maps(G1, G2)
  g, origins = _product_graph_and_edges(G1, G2)
  product_flags = Vector{Vector{OrbifoldFlagWeight{R}}}(undef, n1 * n2)
  vertex_isotropy = Vector{OrbifoldVertexIsotropy}(undef, n1 * n2)
  flag_isotropy = Vector{Vector{OrbifoldFlagIsotropy}}(undef, n1 * n2)
  for v2 in 1:n2, v1 in 1:n1
    V = _product_vertex_index(v1, v2, n1)
    product_flags[V] = vcat(
      [OrbifoldFlagWeight{R}(f1(weight(G1, v1, i)),
        order_of_generic_stabilizer(flags(G1, v1)[i])) for i in 1:val1],
      [OrbifoldFlagWeight{R}(f2(weight(G2, v2, i)),
        order_of_generic_stabilizer(flags(G2, v2)[i])) for i in 1:val2],
    )
    I1, I2 = G1.vertex_isotropy[v1], G2.vertex_isotropy[v2]
    vertex_isotropy[V] = OrbifoldVertexIsotropy(
      vcat(I1.isotropy_group, I2.isotropy_group),
      _block_diagonal(I1.tangent_rep, I2.tangent_rep),
    )
    flag_isotropy[V] = vcat(
      [OrbifoldFlagIsotropy(vcat(F.isotropy_group, I2.isotropy_group),
        _block_diagonal(F.embedding, identity_matrix(ZZ, length(I2.isotropy_group))))
       for F in G1.flag_isotropy[v1]],
      [OrbifoldFlagIsotropy(vcat(I1.isotropy_group, F.isotropy_group),
        _block_diagonal(identity_matrix(ZZ, length(I1.isotropy_group)), F.embedding))
       for F in G2.flag_isotropy[v2]],
    )
  end
  product_labels = _product_labels(G1, G2)
  core = GKMCombinatorialData{R,eltype(product_labels),OrbifoldFlagWeight{R}}(
    g, M, product_labels, product_flags, _product_edge_flags(G1, G2, g, origins),
  )
  con = _product_connection(G1, G2, origins, val1, val2, QQFieldElem)
  return OrbifoldGKMGraph(core, vertex_isotropy, flag_isotropy, con)
end
