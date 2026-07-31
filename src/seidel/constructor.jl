@doc raw"""
    Seidel_space(G::GKMGraph{R}, parameter::AbstractAlgebra.Generic.FreeModuleElem{R}; basePoint::Int=1) where {R} -> GKMGraph

Construct the Seidel space associated to the GKM graph `G` (representing the GKM variety $X$) and the map $\iota:\mathbb{C}^\times\rightarrow T$ given by the element $w\in\mathfrak{t}$.

# Optional arguments:
 - `basePoint::Int64`: This is the vertex of `G` so that in the internal presentation of the curve classes of $S_X$, the curve class of the section of $S_X\rightarrow \mathbb{P}^1$
    associated to the vertex `basePoint` is represented by `(0,...,0,1)`.
    The first entries correspond to curve classes of $X$. The last is the degree of the curve class projected to $\mathbb{P}^1$.

# Examples
```jldoctest Seidel_space
julia> G = projective_space(GKMGraph, 2);

julia> M = lattice(G) # character lattice of G
Free module of rank 3 over ZZ

julia> gens_of_M = gens(M)
3-element Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}:
 (1, 0, 0)
 (0, 1, 0)
 (0, 0, 1)

julia> S = Seidel_space(G, gens_of_M[1])
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, 1)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 0)
Seidel(Birkhoff-Grothendieck) connection for GKM graph with 6 nodes and valency 3

julia> S = Seidel_space(G, gens_of_M[1] + 7*gens_of_M[2])
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, -6)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 7)
Seidel(Birkhoff-Grothendieck) connection for GKM graph with 6 nodes and valency 3

julia> print_curve_classes(S)
[2]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [2]_0: (1, 0), Chern number: 3
[1]_inf -> [1]_0: (0, 1), Chern number: -3
[2]_inf -> [2]_0: (6, 1), Chern number: 15
[2]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [3]_0: (-1, 1), Chern number: -6
[3]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [2]_inf: (1, 0), Chern number: 3

```
Using a different base point does not change the resulting GKM graph but gives a different internal presentation of the curve classes.
```jldoctest Seidel_space
julia> S = Seidel_space(G, gens_of_M[1] + 7*gens_of_M[2]; basePoint=2)
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, -6)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 7)
Seidel(Birkhoff-Grothendieck) connection for GKM graph with 6 nodes and valency 3

julia> print_curve_classes(S)
[2]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [2]_0: (1, 0), Chern number: 3
[1]_inf -> [1]_0: (-6, 1), Chern number: -3
[2]_inf -> [2]_0: (0, 1), Chern number: 15
[2]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [3]_0: (-7, 1), Chern number: -6
[3]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [2]_inf: (1, 0), Chern number: 3
```
"""
function Seidel_space(G::GKMGraph{R}, parameter::AbstractAlgebra.Generic.FreeModuleElem{R}; basePoint::Int=1) where {R}
  @req parent(parameter) === lattice(G) "The weight does not belong to the right character lattice"
  @req is_connected(graph(G)) "The GKM graph should be connected."
  @req 1 <= basePoint <= num_vertices(G) "basePoint is not a vertex of G"
  n, r, d = num_vertices(G), rank_torus(G), valency(G)
  M = free_module(base_ring(lattice(G)), r + 1)
  z = gens(M)[end]
  inclusion = ModuleHomomorphism(lattice(G), M, gens(M)[1:r])
  new_labels = vcat([Vertex("[$(label(G,v))]_0") for v in 1:n], [Vertex("[$(label(G,v))]_inf") for v in 1:n])
  new_flags = [FlagWeight{R}[] for _ in 1:2n]
  for v in 1:n, i in 1:d
    alpha = weight(G, v, i)
    push!(new_flags[v], FlagWeight{R}(inclusion(alpha)))
    pairing = sum(alpha[j] * parameter[j] for j in 1:r)
    push!(new_flags[n+v], FlagWeight{R}(inclusion(alpha) - pairing*z))
  end
  for v in 1:n
    push!(new_flags[v], FlagWeight{R}(z)); push!(new_flags[n+v], FlagWeight{R}(-z))
  end
  new_graph = Graph{Undirected}(2n)
  edge_flags = Dict{Edge,Tuple{Int,Int}}()
  function connect!(u, v, i, j)
    add_edge!(new_graph, u, v)
    stored = Edge(u,v) in edges(new_graph) ? Edge(u,v) : Edge(v,u)
    edge_flags[stored] = src(stored) == u ? (i,j) : (j,i)
  end
  for e in edges(G)
    i, j = core(G).edge_flags[e]
    connect!(src(e), dst(e), i, j); connect!(src(e)+n, dst(e)+n, i, j)
  end
  for v in 1:n; connect!(v, n+v, d+1, d+1); end
  data = GKMCombinatorialData{R,Vertex,FlagWeight{R}}(new_graph, M, new_labels, new_flags, edge_flags)

  old_H2 = GKM_second_homology(G)
  edge_list, edge_to_gen = _edge_list_and_indices(data)
  edge_lattice = free_module(ZZ, length(edge_list))
  H2 = free_module(ZZ, rank(old_H2.H2)+1)
  section = gens(H2)[end]
  embed(beta) = H2(vcat([beta[i] for i in 1:rank(old_H2.H2)], [ZZ(0)]))
  shifts = Dict{Int,CurveClass}(1 => zero(old_H2.H2))
  while length(shifts) < n
    for v in collect(keys(shifts)), u in all_neighbors(graph(G), v)
      haskey(shifts,u) && continue
      e = Edge(v,u); alpha = weight(G,e)
      shifts[u] = shifts[v] - sum(alpha[i]*parameter[i] for i in 1:r)*curve_class(G,e)
    end
  end
  images = Vector{typeof(zero(H2))}(undef, length(edge_list))
  for (i,e) in enumerate(edge_list)
    if src(e)<=n && dst(e)<=n
      images[i] = embed(curve_class(G,e))
    elseif src(e)>n && dst(e)>n
      images[i] = embed(curve_class(G,Edge(src(e)-n,dst(e)-n)))
    else
      v = min(src(e),dst(e)); images[i] = embed(shifts[v]-shifts[basePoint]) + section
    end
  end
  quotient = ModuleHomomorphism(edge_lattice, H2, images)
  dual_cone, ray_sum, chern = _finish_GKM_H2(data, edge_lattice, H2, quotient, edge_list)
  seidel_H2 = GKM_H2(edge_lattice, H2, edge_to_gen, quotient, dual_cone, ray_sum, chern)
  # Natural connection: copy the connection on G to both fibres, keep the
  # section flag parallel along fibre edges, and identify corresponding flags
  # along each section edge.
  seidel_transport = Dict{Edge,Vector{Int}}()
  for unoriented_edge in edges(G), e in (unoriented_edge, reverse(unoriented_edge))
    fibre_transport = vcat(transport(G)[e], d + 1)
    seidel_transport[e] = fibre_transport
    seidel_transport[Edge(src(e) + n, dst(e) + n)] = copy(fibre_transport)
  end
  for v in 1:n
    section_edge = Edge(v, n + v)
    seidel_transport[section_edge] = collect(1:d+1)
    seidel_transport[reverse(section_edge)] = collect(1:d+1)
  end
  seidel_coefficients = Dict{Edge,Vector{R}}()
  for (e, image) in seidel_transport
    seidel_coefficients[e] = _coefficients_for_transport(data, e, image, R)
  end
  con = Connection{R}(
    seidel_transport, seidel_coefficients, "Seidel($(connection_type(G)))",
  )
  return GKMGraph{R,Vertex,FlagWeight{R}}(data, con, create_cohomology(r+1,2n), seidel_H2, nothing)
end

const _seidel_section_counts = IdDict{Any,Any}()
function _SeidelSectionCount(SG::GKMGraph)
  @req iseven(num_vertices(SG)) "SG is not a Seidel space!"
  get!(_seidel_section_counts, SG) do
    H2 = GKM_second_homology(SG).H2; Z = free_module(ZZ,1)
    ModuleHomomorphism(H2, Z, [i==rank(H2) ? gens(Z)[1] : zero(Z) for i in 1:rank(H2)])
  end
end

function _effectiveSectionClassesWithChernNumber(SG::GKMGraph, chernNumber::ZZRingElem)
  H2 = GKM_second_homology(SG); sec = _SeidelSectionCount(SG)
  target, _, _ = direct_sum([codomain(H2.chern), codomain(sec)])
  q = ModuleHomomorphism(H2.H2, target, hcat(matrix(H2.chern),matrix(sec)))
  success, e0 = has_preimage_with_preimage(q, chernNumber*gens(target)[1]+gens(target)[2])
  success || return Vector{}()
  K, k = kernel(q); mk = transpose(matrix(k)); dual_rays = rays(H2.dual_cone)
  ray_matrix = QQMatrix(length(dual_rays),rank(H2.H2))
  for i in eachindex(dual_rays), j in 1:rank(H2.H2); ray_matrix[i,j]=dual_rays[i][j]; end
  P = polyhedron(-ray_matrix*mk, ray_matrix*[e0[i] for i in 1:rank(H2.H2)])
  return (e0+k(K([p[i] for i in 1:rank(K)])) for p in lattice_points(P))
end
