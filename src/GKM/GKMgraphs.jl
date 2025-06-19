@doc raw"""
    gkm_graph(g, labels, M, w; check=true, checkLabels=true) -> GKM_graph

Create a GKM graph from the given data.

# Arguments
- `g::Graph`: An unoriented OSCAR graph.
- `labels::Vector{String}`: A vector of strings, used to denote the vertices.
- `M::AbstractAlgebra.Generic.FreeModule{R}`: A OSCAR free module over ``\mathbb{Z}`` or ``\mathbb{Q}``, it denotes the character group. `R` is either `ZZRingElem` or `QQFieldElem`.
- `w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}`: The axial function. Note that it is enough to specify the weight of each edge in one orientation here. The opposite oriented edge will automatically be given minus that weight.
- `check::Bool=true`: Check if the data inserted are consistent.
- `checkLabels::Bool=true`: Check that the labels don't contain the characters `<`, `[`, `]`, which are reserved for the output of special constructions like blowups and projective bundles.

# Example
Let us construct the GKM graph of the projective line. First of all, we create a graph with two vertices, and one edge.
```jldoctest first_GKM_graph
julia> g = Graph{Undirected}(2)
Undirected graph with 2 nodes and no edges

julia> add_edge!(g, 1, 2);
```
Let us define our array of labels.

```jldoctest first_GKM_graph
julia> labels = ["a", "b"];
```
Now, we create the character group. We take a free module of rank 2 over the integers.

```jldoctest first_GKM_graph
julia> M = free_module(ZZ, 2)
Free module of rank 2 over ZZ
```
We create the axial function. It is a dictionary from the set of edges to the character group. This time we have only one edge.

```jldoctest first_GKM_graph
julia> e = first(edges(g));

julia> w = Dict(e => gens(M)[1] - gens(M)[2])
Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}} with 1 entry:
  Edge(2, 1) => (1, -1)
```
Finally, we create the GKM graph.
```jldoctest first_GKM_graph
julia> gkm_graph(g, labels, M, w)
GKM graph with 2 nodes, valency 1 and axial function:
b -> a => (1, -1)
```

!!! warning
    1. Do not change the number of verices after this.
    2. Don't modify the underlying OSCAR graph directly after this. Use the functions of this package instead.
    3. All edges should be added immediately after calling this function and not changed afterwards.

!!! note
    After you have added all edges using `add_edge!`, you may use `initialize!` to calculate the GKM connection (if it is unique) and the curve classes. If you don't do this, those data will be calculated whenever required for the first time.

"""
function gkm_graph(
  g::Graph,
  labels::Vector{String},
  M::AbstractAlgebra.Generic.FreeModule{R}, # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}; 
  check::Bool=true,
  checkLabels::Bool=true
) where R <: Weight
  # construct the GKM_graph
  if check
    @req n_vertices(g) >= 1 "GKM graph needs at least one vertex"
    @req length(labels) == n_vertices(g) "The number of labels does not match the number of fixed points"
    @req all(e -> parent(w[e]) === M, edges(g)) "Character group mismatch"
    @req Set(edges(g)) == keys(w) "The axial function is not well defined"
    @req all(v -> length(all_neighbors(g, 1)) == length(all_neighbors(g, v)), 2:n_vertices(g)) "The valency is not the same for all vertices"
    @req length(unique(labels)) == length(labels) "Labels must be unique"
  end

  if checkLabels
    # reserve characters <,[,] for vertex labels of blowups and Seidel space
    @req all(v -> !contains(labels[v], ">") && !contains(labels[v], "[") && !contains(labels[v], "]"), 1:n_vertices(g)) "Characters >,[,] are forbidden for vertex labels"
  end

  for e in edges(g)
    w[reverse(e)] = -w[e]
  end

  QH_structure_consts = Dict{CurveClass, Array{Any, 3}}()

  coeffRing, _ = polynomial_ring(QQ, ["t$i" for i in 1:rank(M)])
  coeffRingLocalized = fraction_field(coeffRing)
  cohomRing = free_module(coeffRing, n_vertices(g))
  cohomRingLocalized = free_module(coeffRingLocalized, n_vertices(g))
  edgeWeightClasses = Dict{Edge, QQMPolyRingElem}()
  pointEulerClasses::Vector{Union{Nothing, QQMPolyRingElem}} = [nothing for _ in 1:n_vertices(g)]#vcat(Union{Nothing, QQMPolyRingElem}[], repeat([nothing], n_vertices(G.g)))

  coh_ring = GKM_cohomology_ring(coeffRing, coeffRingLocalized, cohomRing, cohomRingLocalized, edgeWeightClasses, pointEulerClasses)

  gkm = GKM_graph(g, labels, M, w, coh_ring, nothing, nothing, QH_structure_consts, false)

  return gkm
end


@doc raw"""
    empty_gkm_graph(n::Int64, r::Int64, labels::Vector{String}) -> GKM_graph

Return the GKM graph with `n` fixed points, no edges, torus rank `r` and vertices labelled by `labels`.

# Example

```jldoctest empty_GKM_graph
julia> G = empty_gkm_graph(2, 2, ["a", "b"])
GKM graph with 2 nodes, valency 0 and axial function:

```
"""
function empty_gkm_graph(n::Int64, r::Int64, labels::Vector{String})

  return gkm_graph(Graph{Undirected}(n), labels, free_module(ZZ, r), Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}())
end

@doc raw"""
    add_edge!(G::GKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R <: Weight

Add an edge to `G` from the vertex labelled `s` to the vertex labelled `d`, the axial function takes value `weight` in that edge.

Let us construct the same example of `gkm_graph`, that is the GKM graph of the projective space. 

# Example

```jldoctest empty_GKM_graph
julia> G = empty_gkm_graph(2, 2, ["a", "b"])
GKM graph with 2 nodes, valency 0 and axial function:

julia> wei = gens(G.M)[1] - gens(G.M)[2]
(1, -1)

julia> add_edge!(G, "b", "a", wei);

julia> G
GKM graph with 2 nodes, valency 1 and axial function:
b -> a => (1, -1)

```

!!! warning    
    Add all edges immediately after creation. Any curve class or cohomology functionality should only be used after all edges have been added. The same holds for `initialize!`.
"""
function add_edge!(G::GKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R <: Weight

  @req (s in G.labels) "Source label not found"
  @req (d in G.labels) "Destination label not found"

  sd = indexin([s, d], G.labels)

  add_edge!(G, sd[1], sd[2], weight)

end

@doc raw"""
    add_edge!(G::GKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R <: Weight

Same as before, but using the number of the vertex instead of the label.
"""
function add_edge!(G::GKM_graph, s::Int64, d::Int64, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R <: Weight
  
  @req (s in 1:n_vertices(G.g)) "Source not found"
  @req (d in 1:n_vertices(G.g)) "Destination not found"
  @req parent(weight) === G.M "The group of characters is not correct"

  Oscar.add_edge!(G.g, s, d)
  G.w[Edge(s, d)] = weight
  G.w[Edge(d, s)] = -weight

end

@doc raw"""
    valency(G::GKM_graph) -> Int64

Return the valency of `G`, i.e. the degree of each vertex.
# Example:
The valency of the GKM graph of $\mathbb{P}^3$ is 3, since all of the fixed points $[1:0:0:0], \dots, [0:0:0:1]$ are connected to each other
via some $T$-invariant $\mathbb{P}^1$'s. For example, $[1:0:0:0]$ and $[0:1:0:0]$ are connected by $\{[x:y:0:0] : x,y\in\mathbb{C}\}$.
```jldoctest valency
julia> valency(projective_space(GKM_graph, 3))
3
julia> valency(grassmannian(GKM_graph, 2, 4)) # The Grassmannian of 2-planes in C^4
4
julia> valency(flag_variety(GKM_graph, [1, 1, 1, 1])) # The variety of full flags in C^4
6
```
"""
function valency(G::GKM_graph)
  return length(all_neighbors(G.g, 1))
end


@doc raw"""
    rank_torus(G::GKM_graph) -> Int64

Return the rank of the torus acting on `G`. That is, the rank of the character group.

# Examples
By default, the torus acting on $\mathbb{P}^n$ is $(\mathbb{C}^\times)^{n+1}$, acting by rescaling the homogeneous coordinates.
```jldoctest rank_torus
julia> P3 = projective_space(GKM_graph, 3);

julia> rank_torus(P3)
4
```
Taking products adds the rank:
```jldoctest rank_torus
julia> H6 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 6));

julia> rank_torus(H6)
4
julia> rank_torus(H6 * P3)
8
```
"""
function rank_torus(G::GKM_graph)
  return rank(G.M)
end

@doc raw"""
    is2_indep(G::GKM_graph) -> Bool

Return `true` if `G` is 2-independent, i.e. the weights of every two edges at a vertex are linearly independent.
"""
function is2_indep(G::GKM_graph)
  return _indep(G, 2)
end

@doc raw"""
    is3_indep(G::GKM_graph) -> Bool

Return `true` if `G` is 3-independent, i.e. the weights of every three edges at a vertex are linearly independent.
# Example
The weights of $\mathbb{P}^3$ at the fixed point $[1:0:0:0]$ are $\{t_i-t_0:i\in\{1, 2, 3\}\}$, which are linearly independent over $\mathbb{C}$.
```jldoctest is3_indep
julia> is3_indep(projective_space(GKM_graph, 3))
true
```
The variety of complete flags in $\mathbb{C}^3$ is an example of a GKM graph that is not 3-independent:
```jldoctest is3_indep
julia> G = flag_variety(GKM_graph, [1, 1, 1])
GKM graph with 6 nodes, valency 3 and axial function:
13 -> 12 => (0, -1, 1)
21 -> 12 => (-1, 1, 0)
23 -> 13 => (-1, 1, 0)
23 -> 21 => (-1, 0, 1)
31 -> 13 => (-1, 0, 1)
31 -> 21 => (0, -1, 1)
32 -> 12 => (-1, 0, 1)
32 -> 23 => (0, -1, 1)
32 -> 31 => (-1, 1, 0)

julia> is3_indep(G)
false
```
!!! warning
    This function throws an error if the valency of `G` is less than 3, since in this case it is not possible to pick three different edges ta a vertex.
"""
function is3_indep(G::GKM_graph)
  return _indep(G, 3)
end

function _indep(G::GKM_graph, k::Int64)
  
  @req valency(G) >= k "valency is too low"

  for v in 1:n_vertices(G.g)
    for tup in Iterators.product([all_neighbors(G.g, v) for _ in 1:k]...)
      any(i-> tup[i-1] >= tup[i], 2:k) && continue
      
      if rank(matrix([G.w[Edge(v, tup[i])] for i in 1:k])) < k
        return false
      end

    end
  end
  
  return true
end

function Base.show(io::IO, G::GKM_graph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM graph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM graph with $(n_vertices(G.g)) nodes and valency $(valency(G))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::GKM_graph)

  print(io, "GKM graph with $(n_vertices(G.g)) nodes, valency $(valency(G)) and axial function:")
  for e in edges(G.g)
    print(io, "\n$(G.labels[src(e)]) -> $(G.labels[dst(e)]) => $(G.w[e])")
  end
end

"""
    initialize!(gkm::GKM_graph; connection::Bool=true, curveClasses::Bool=true)

You may optionally call this function as soon as all edges have been added to the GKM graph to calculate the GKM connection
(if unique) and the curve classes of the gkm graph.
This will set the fields `gkm.connection` and `gkm.H2` that are initially `nothing`.
If you don't call this, these fields will be initialized later if possible, which might take some time at unexpected moments
(especially for curve classes).
If any of those fields are already set, this will not overwrite them.
"""
function initialize!(gkm::GKM_graph; connection::Bool=true, curveclasses::Bool=true)

  if connection
    get_connection(gkm)
  end
  if curveclasses
    GKM_second_homology(gkm)
  end
end


function _common_weight_denominator(G::GKM_graph)::ZZRingElem
  
  if isa(_get_weight_type(G), ZZRingElem)
    return ZZ(1)
  elseif isa(_get_weight_type(G), QQFieldElem)
    res::ZZRingElem = ZZ(1)
    rk = rank_torus(G)
    for e in edges(G.g)
      for i in 1:rk
        res = lcm(res, denominator(G.w[e][i]))
      end
    end
    return res
  end
  @req false "Ony ZZRingElem and QQFieldElem are supported as weight types so far."
end

@doc raw"""
    isvalid(gkm::GKM_graph; printDiagnostics::Bool=true) -> Bool
Return true if the GKM graph is valid. This means:
  1. Every vertex has the same degree
  2. The weights are defined for every edge and every reverse of every edge
  3. The weights belong to the weight lattice
  4. The weights of an edge and its reverse sum to zero
  5. There are the right number of vertex labels
  6. If the valency is at least two, the weights of the graph are 2-independent.
  7. Vertex labels must be unique
  8. The equivariant cohomology ring has rank = number of vertices of graph
  9. The coefficient ring of the equivariant cohomology ring has number of generators = torus rank.

# Examples
The standard constructions always produce valid GKM graphs, e.g. the complex projective space $\mathbb{P}^3$:
```jldoctest isvalid_GKM_graph
julia> isvalid(projective_space(GKM_graph, 3))
true
```
On the other hand, here is an example showing why one should never modify the underlying OSCAR graph of a GKM graph directly:
```jldoctest isvalid_GKM_graph
julia> G = empty_gkm_graph(3, 1, ["v1", "v2", "v3"])
GKM graph with 3 nodes, valency 0 and axial function:

julia> add_edge!(G.g, 1, 2)
true

julia> isvalid(G)
The valency is not the same for all vertices
false

julia> add_edge!(G.g, 1, 3)
true

julia> add_edge!(G.g, 2, 3)
true

julia> isvalid(G)
Weight of Edge(2, 1) is missing.
false
```
Instead, one should add all edges with `add_edge!(G, "v1", "v2", weight)` (see above).
"""
function isvalid(gkm::GKM_graph; printDiagnostics::Bool=true)::Bool

  if !all(v -> length(all_neighbors(gkm.g, 1)) == length(all_neighbors(gkm.g, v)), 2:n_vertices(gkm.g))
    printDiagnostics && println("The valency is not the same for all vertices")
    return false
  end

  for e in edges(gkm.g)
    if !haskey(gkm.w, e)
      printDiagnostics && println("Weight of $e is missing.")
      return false
    elseif !haskey(gkm.w, reverse(e))
      printDiagnostics && println("Weight of $(reverse(e)) is missing.")
      return false
    elseif parent(gkm.w[e]) != gkm.M
      printDiagnostics && println("Weight of $e doesn't belong to $(gkm.M).")
      return false
    elseif parent(gkm.w[reverse(e)]) != gkm.M
      printDiagnostics && println("Weight of $(reverse(e)) doesn't belong to $(gkm.M).")
      return false
    elseif !(gkm.w[e] == -gkm.w[reverse(e)])
      printDiagnostics && println("Weights of $e and $(reverse(e)) don't sum to zero.")
      return false
    end
  end

  if length(gkm.labels) != n_vertices(gkm.g)
    printDiagnostics && println("Not the right number of labels")
    return false
  elseif (valency(gkm) > 1 && !is2_indep(gkm))
    printDiagnostics && println("GKM graph is not 2-independent.")
    return false
  end

  if length(unique(gkm.labels)) != length(gkm.labels)
    printDiagnostics && println("Labels are not unique.")
    return false
  end

  if (valency(gkm) > 2 && !is3_indep(gkm))
    printDiagnostics && println("GKM graph is valid but not 3-independent, so connections may not be unique.")
  end

  if length(gens(gkm.equivariantCohomology.coeffRing)) != rank_torus(gkm)
    printDiagnostics && println("Coefficient ring of equivariant cohomology has wrong rank.")
    return false
  end

  if length(gens(gkm.equivariantCohomology.cohomRing)) != n_vertices(gkm.g)
    printDiagnostics && println("Equivariant cohomology ring should have rank = number of vertices.")
    return false
  end

  return true
end

function edgeFromLabels(G::GKM_graph, s::String, d::String)::Edge
  @req s in G.labels "source $s is not a vertex in $G."
  @req d in G.labels "destination $d is not a vertex in $G."

  sd = indexin([s, d], G.labels)
  return Edge(sd[1], sd[2])
end
