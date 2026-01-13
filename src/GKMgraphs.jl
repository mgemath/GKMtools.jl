@doc raw"""
    gkm_graph(g, labels, M, w; check=true, checkLabels=true) -> AbstractGKM_graph
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
) where R <: GKM_weight_type
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

  # Build flag-based data structures from edge weights
  nv = n_vertices(g)
  weights_at_vertex = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}}(undef, nv)
  flag_to_edge = Vector{Vector{Union{Nothing, Edge}}}(undef, nv)
  edge_to_flag_index = Dict{Edge, Int64}()

  for v in 1:nv
    neighbors = all_neighbors(g, v)
    deg = length(neighbors)
    weights_at_vertex[v] = Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}(undef, deg)
    flag_to_edge[v] = Vector{Union{Nothing, Edge}}(undef, deg)

    for (i, nbr) in enumerate(neighbors)
      e = Edge(v, nbr)
      weights_at_vertex[v][i] = w[e]
      flag_to_edge[v][i] = e
      edge_to_flag_index[e] = i
    end
  end

  GW_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()

  gkm = AbstractGKM_graph(g, labels, M, weights_at_vertex, edge_to_flag_index, flag_to_edge, w, nothing, nothing, nothing, GW_structure_consts, false)
  gkm.equivariantCohomology = _equivariant_cohomology_ring(gkm)

  return gkm
end

@doc raw"""
    flags_only_gkm_graph(labels::Vector{String}, M::AbstractAlgebra.Generic.FreeModule{R}, w::Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}}; check::Bool=true, checkLabels::Bool=true) -> AbstractGKM_graph

Same as `gkm_graph`, but constructs an empty GKM graph without any edges and only standalone flags.
The input `w[v][i]` is the weight of the `i`-th flag at vertex `v`.
Edges may be added later using `connect_flags!`.

# Example
```jldoctest
julia> M = free_module(ZZ, 2);

julia> g1, g2 = gens(M);

julia> G = flags_only_gkm_graph(["v1", "v2"], M, [[g1, g2], [-g1, g1+g2]])
GKM graph with 2 nodes, valency 2 and axial function:
Standalone flags:
v1.1 => (1, 0)
v1.2 => (0, 1)
v2.1 => (-1, 0)
v2.2 => (1, 1)

julia> connect_flags!(G, 1, 2, 1, 1);

julia> G
GKM graph with 2 nodes, valency 2 and axial function:
v2 -> v1 => (-1, 0)
Standalone flags:
v1.2 => (0, 1)
v2.2 => (1, 1)
```
"""
function flags_only_gkm_graph(labels::Vector{String},
  M::AbstractAlgebra.Generic.FreeModule{R}, # character group
  w::Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}};
  check::Bool=true,
  checkLabels::Bool=true) where R <: GKM_weight_type

  nv = length(w)
  deg = length(w[1])

  if check
    @req nv >= 1 "GKM graph needs at least one vertex"
    @req length(labels) == nv "The number of labels does not match the number of fixed points"
    @req all(w_at_flag -> all(wt -> parent(wt) === M, w_at_flag), w) "Character group mismatch"
    @req all(v -> length(w[1]) == length(w[v]), 2:nv) "The valency is not the same for all vertices"
    @req length(unique(labels)) == length(labels) "Labels must be unique"
  end
  if checkLabels
    # reserve characters <,[,] for vertex labels of blowups and Seidel space
    @req all(v -> !contains(labels[v], ">") && !contains(labels[v], "[") && !contains(labels[v], "]"), 1:nv) "Characters >,[,] are forbidden for vertex labels"
  end

  weights_at_vertex = w
  flag_to_edge = Vector{Vector{Union{Nothing, Edge}}}(undef, nv)
  edge_to_flag_index = Dict{Edge, Int64}()

  for v in 1:nv
    flag_to_edge[v] = Vector{Union{Nothing, Edge}}(nothing, deg)
  end

  GW_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()
  g = Graph{Undirected}(nv)
  w_old = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}()

  gkm = AbstractGKM_graph(g, labels, M, weights_at_vertex, edge_to_flag_index, flag_to_edge, w_old, nothing, nothing, nothing, GW_structure_consts, false)
  gkm.equivariantCohomology = _equivariant_cohomology_ring(gkm)

  return gkm
end

@doc raw"""
    connect_flags!(G::AbstractGKM_graph, v1::Int64, v2::Int64, f1::Int64, f2::Int64) -> Edge

Connect flag `f1` at vertex `v1` and flag `f2` at vertex `v2` to form a new edge of `G`, 
and return the new edge.
This requires that vertices `v1` and `v2` are not yet connected by an edge,
and that the weights of the two flags sum to zero.

!!! warning    
    Add all edges immediately after creation. Any curve class or cohomology functionality should only be used after all edges have been added. The same holds for `initialize!`.

An example of this function is provided in `flags_only_gkm_graph` above.
"""
function connect_flags!(G::AbstractGKM_graph, v1::Int64, v2::Int64, f1::Int64, f2::Int64)
  nv = n_vertices(G.g)
  val = valency(G)
  @req (1 <= v1) && (v1 <= nv) "Vertex v1=$v1 out of bounds."
  @req (1 <= v2) && (v2 <= nv) "Vertex v2=$v2 out of bounds."
  @req (1 <= f1) && (f1 <= val) "Flag index f1=$f1 out of bounds."
  @req (1 <= f2) && (f2 <= val) "Flag index f2=$f2 out of bounds."
  @req isnothing(G.flag_to_edge[v1][f1]) "Flag $f1 at vertex $v1 already belongs to an edge."
  @req isnothing(G.flag_to_edge[v2][f2]) "Flag $f2 at vertex $v2 already belongs to an edge."
  @req iszero(G.weights_at_vertex[v1][f1] + G.weights_at_vertex[v2][f2]) "Flag weights must sum to zero."
  @req !has_edge(G.g, v1, v2) "Edge($v1, $v2) already exists in GKM graph."

  e = Edge(v1, v2)
  w = G.weights_at_vertex[v1][f1]
  add_edge!(G.g, v1, v2)
  G.edge_to_flag_index[e] = f1
  G.edge_to_flag_index[reverse(e)] = f2
  G.flag_to_edge[v1][f1] = e
  G.flag_to_edge[v2][f2] = reverse(e)
  G.w[e] = w
  G.w[reverse(e)] = -w

  return e
end

@doc raw"""
    empty_gkm_graph(n::Int64, r::Int64, labels::Vector{String}) -> AbstractGKM_graph

It returns the GKM graph with `n` fixed points, no edges, torus rank `r` and vertices labelled by `labels`.

```jldoctest empty_GKM_graph
julia> G = empty_gkm_graph(2, 2, ["a", "b"])
GKM graph with 2 nodes, valency 0 and axial function:

```
"""
function empty_gkm_graph(n::Int64, r::Int64, labels::Vector{String})

  return gkm_graph(Graph{Undirected}(n), labels, free_module(ZZ, r), Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}())
end

@doc raw"""
    add_edge!(G::AbstractGKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) -> Tuple{Int64, Int64} where R<:GKM_weight_type

Add an edge to `G` from the vertex labelled `s` to the vertex labelled `d`, the axial function takes value `weight` in that edge.
The returned tuple `(i,j)` gives the indices of the two flags of the newly created edge at the vertex `s` and `d`, respectively.

Let us construct the same example of `gkm_graph`, that is the GKM graph of the projective space. 

```jldoctest empty_GKM_graph
julia> G = empty_gkm_graph(2, 2, ["a", "b"])
GKM graph with 2 nodes, valency 0 and axial function:

julia> wei = gens(G.M)[1] - gens(G.M)[2]
(1, -1)

julia> add_edge!(G, "b", "a", wei)
(1, 1)

julia> G
GKM graph with 2 nodes, valency 1 and axial function:
b -> a => (1, -1)

```

!!! warning    
    Add all edges immediately after creation. Any curve class or cohomology functionality should only be used after all edges have been added. The same holds for `initialize!`.
"""
function add_edge!(G::AbstractGKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

  @req (s in G.labels) "Source label not found"
  @req (d in G.labels) "Destination label not found"

  sd = indexin([s, d], G.labels)

  add_edge!(G, sd[1], sd[2], weight)

end

@doc raw"""
    add_edge!(G::AbstractGKM_graph, s::String, d::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

Same as before, but using the number of the vertex instead of the label.
"""
function add_edge!(G::AbstractGKM_graph, s::Int64, d::Int64, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

  @req (s in 1:n_vertices(G.g)) "Source $s not found"
  @req (d in 1:n_vertices(G.g)) "Destination $d not found"
  @req parent(weight) === G.M "The group of characters is not correct"

  Oscar.add_edge!(G.g, s, d)

  # Update backward-compatible w dict
  G.w[Edge(s, d)] = weight
  G.w[Edge(d, s)] = -weight

  # Update flag-based structures
  # Add flag at vertex s for edge (s,d)
  push!(G.weights_at_vertex[s], weight)
  push!(G.flag_to_edge[s], Edge(s, d))
  G.edge_to_flag_index[Edge(s, d)] = length(G.weights_at_vertex[s])

  # Add flag at vertex d for edge (d,s)
  push!(G.weights_at_vertex[d], -weight)
  push!(G.flag_to_edge[d], Edge(d, s))
  G.edge_to_flag_index[Edge(d, s)] = length(G.weights_at_vertex[d])
  return (length(G.weights_at_vertex[s]), length(G.weights_at_vertex[d]))
end

@doc raw"""
    add_standalone_flag!(G::AbstractGKM_graph, v::Int64, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) -> Int64 where R<:GKM_weight_type

Add a standalone flag (not associated with an edge) at vertex `v` with the given `weight`
and return the index of the newly created flag at `v`.

# Example
```jldoctest add_standalone_flag
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> g1, g2, g3 = gens(G.M);

julia> add_standalone_flag!(G, 1, g1)
3

julia> add_standalone_flag!(G, 2, g2);

julia> add_standalone_flag!(G, 3, g3);

julia> valency(G)
3

julia> G
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Standalone flags:
1.3 => (1, 0, 0)
2.3 => (0, 1, 0)
3.3 => (0, 0, 1)
```
"""
function add_standalone_flag!(G::AbstractGKM_graph, v::Int64, weight::AbstractAlgebra.Generic.FreeModuleElem{R})::Int64 where R<:GKM_weight_type

  @req (v in 1:n_vertices(G.g)) "Vertex $v not found"
  @req parent(weight) === G.M "The group of characters is not correct"

  # Add the standalone flag
  push!(G.weights_at_vertex[v], weight)
  push!(G.flag_to_edge[v], nothing)  # No associated edge

  return length(G.weights_at_vertex[v])
end

@doc raw"""
    add_standalone_flag!(G::AbstractGKM_graph, v::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

Same as above, but using vertex label instead of index.
"""
function add_standalone_flag!(G::AbstractGKM_graph, v::String, weight::AbstractAlgebra.Generic.FreeModuleElem{R}) where R<:GKM_weight_type

  @req (v in G.labels) "Vertex label $v not found"

  v_idx = findfirst(==(v), G.labels)
  add_standalone_flag!(G, v_idx, weight)
end

@doc raw"""
    valency(G::AbstractGKM_graph) -> Int64

Return the valency of `G`, i.e. the number of flags at each vertex.

!!! warning
    This function does not check if `G` is a valid GKM graph (use `isvalid` to check this).
    In particular, it does not check if every vertex has the same degree.
    The returned value is the degree of vertex `1`.

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
function valency(G::AbstractGKM_graph)
  return length(G.weights_at_vertex[1])
end

@doc raw"""
    is_compact(G::AbstractGKM_graph) -> Bool

Return `true` if `G` is compact, i.e. all flags at all vertices are associated with edges (no standalone flags).

# Example
```jldoctest is_compact
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> is_compact(G)
true

julia> M = G.M;

julia> add_standalone_flag!(G, 1, gens(M)[1]);

julia> add_standalone_flag!(G, 2, gens(M)[2]);

julia> add_standalone_flag!(G, 3, gens(M)[3]);

julia> is_compact(G)
false
```
"""
function is_compact(G::AbstractGKM_graph)
  for v in 1:n_vertices(G.g)
    for flag_edge in G.flag_to_edge[v]
      if isnothing(flag_edge)
        return false
      end
    end
  end
  return true
end


@doc raw"""
    rank_torus(G::AbstractGKM_graph) -> Int64

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
function rank_torus(G::AbstractGKM_graph)
  return rank(G.M)
end

@doc raw"""
    is2_indep(G::AbstractGKM_graph) -> Bool

Return `true` if `G` is 2-independent, i.e. the weights of every two flags at a vertex are linearly independent.
"""
function is2_indep(G::AbstractGKM_graph)
  return _indep(G, 2)
end

@doc raw"""
    is3_indep(G::AbstractGKM_graph) -> Bool

Return `true` if `G` is 3-independent, i.e. the weights of every three flags at a vertex are linearly independent.
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
    This function throws an error if the valency of `G` is less than 3, since in this case it is not possible to pick three different flags at a vertex.
"""
function is3_indep(G::AbstractGKM_graph)
  return _indep(G, 3)
end

function _indep(G::AbstractGKM_graph, k::Int64)

  @req valency(G) >= k "valency is too low"

  val = valency(G)

  for v in 1:n_vertices(G.g)
    # Check all k-tuples of distinct flag indices at vertex v
    for tup in Iterators.product([1:val for _ in 1:k]...)
      # Skip if not strictly increasing (to avoid checking same set multiple times)
      any(i -> tup[i-1] >= tup[i], 2:k) && continue

      # Get the weights of the k flags
      weights = [G.weights_at_vertex[v][tup[i]] for i in 1:k]

      if rank(matrix(weights)) < k
        return false
      end
    end
  end

  return true
end

"""
    initialize!(gkm::AbstractGKM_graph; connection::Bool=true, curveClasses::Bool=true)

You may optionally call this function as soon as all edges have been added to the GKM graph to calculate the GKM connection
(if unique) and the curve classes of the gkm graph.
This will set the fields `gkm.connection` and `gkm.curveClasses` that are initially `nothing`.
If you don't call this, these fields will be initialized later if possible, which might take some time at unexpected moments
(especially for curveClasses).
If any of those fields are already set, this will not overwrite them.
"""
function initialize!(gkm::AbstractGKM_graph; connection::Bool=true, curveClasses::Bool=true)

  if connection
    get_connection(gkm)
  end
  if curveClasses
    GKM_second_homology(gkm)
  end
end


function _common_weight_denominator(G::AbstractGKM_graph)::ZZRingElem
  if G.weightType <: ZZRingElem
    return ZZ(1)
  elseif G.weightType <: QQFieldElem
    res::ZZRingElem = ZZ(1)
    rk = rank_torus(G)
    for v in 1:n_vertices(G.g)
      for weight in G.weights_at_vertex[v]
        for i in 1:rk
          res = lcm(res, denominator(weight[i]))
        end
      end
    end
    return res
  end
  @req false "Ony ZZRingElem and QQFieldElem are supported as weight types so far."
end

function Base.show(io::IO, G::AbstractGKM_graph)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM graph")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM graph with $(n_vertices(G.g)) nodes and valency $(valency(G))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", G::AbstractGKM_graph)

  print(io, "GKM graph with $(n_vertices(G.g)) nodes, valency $(valency(G)) and axial function:")
  for e in edges(G.g)
    print(io, "\n$(G.labels[src(e)]) -> $(G.labels[dst(e)]) => $(G.w[e])")
  end
  # print standalone flags if any:
  is_compact(G) && return
  print(io, "\nStandalone flags:")
  for v in 1:n_vertices(G.g)
    for (i, w) in enumerate(G.weights_at_vertex[v])
      !isnothing(G.flag_to_edge[v][i]) && continue
      print(io, "\n$(G.labels[v]).$i => $w")
    end
  end
end



@doc raw"""
    isvalid(gkm::AbstractGKM_graph; printDiagnostics::Bool=true) -> Bool
Return true if the GKM graph is valid. This means:
  1. Every vertex has the same degree (i.e. number of flags).
  2. All flag weights are defined and belong to the weight lattice.
  3. For edges: the flag weight at each end of an edge sums to zero.
  4. Edge-to-flag mappings are consistent.
  5. There are the right number of vertex labels.
  6. If the valency is at least two, the flag weights are 2-independent.
  7. Vertex labels must be unique.
  8. The equivariant cohomology ring has rank = number of vertices of graph.
  9. The coefficient ring of the equivariant cohomology ring has number of generators = torus rank.

# Examples
The standard constructions always produce valid GKM graphs, e.g. the complex projective space $\mathbb{P}^3$:
```jldoctest isvalid_GKM_graph
julia> isvalid(projective_space(GKM_graph, 3))
true
```
"""
function isvalid(gkm::AbstractGKM_graph; printDiagnostics::Bool=true)::Bool

  nv = n_vertices(gkm.g)
  val = valency(gkm)

  # Check all vertices have same number of flags
  if !all(v -> length(gkm.weights_at_vertex[v]) == val, 1:nv)
    printDiagnostics && println("The number of flags is not the same for all vertices")
    return false
  end

  # Check all flag weights are defined and belong to weight lattice
  for v in 1:nv
    if length(gkm.flag_to_edge[v]) != val
      printDiagnostics && println("Vertex $v has inconsistent flag_to_edge length")
      return false
    end

    for i in 1:val
      w = gkm.weights_at_vertex[v][i]
      if parent(w) != gkm.M
        printDiagnostics && println("Weight of flag $i at vertex $v doesn't belong to $(gkm.M).")
        return false
      end
    end
  end

  # Check edge-to-flag mappings and edge weight consistency
  for e_unoriented in edges(gkm.g)
    for e in [e_unoriented, reverse(e_unoriented)]
      if !haskey(gkm.edge_to_flag_index, e)
        printDiagnostics && println("Edge $e missing from edge_to_flag_index.")
        return false
      end

      i_e = gkm.edge_to_flag_index[e]
      i_rev = gkm.edge_to_flag_index[reverse(e)]

      if i_e < 1 || i_e > val
        printDiagnostics && println("Edge $e has invalid flag index $i_e")
        return false
      end

      # Check that flag_to_edge is consistent
      if gkm.flag_to_edge[src(e)][i_e] != e
        printDiagnostics && println("flag_to_edge inconsistency for edge $e")
        return false
      end

      # Check edge weights sum to zero
      w_e = gkm.weights_at_vertex[src(e)][i_e]
      w_rev = gkm.weights_at_vertex[dst(e)][i_rev]
      if !(w_e == -w_rev)
        printDiagnostics && println("Weights of $e and $(reverse(e)) don't sum to zero.")
        return false
      end

      # Check backward compatibility dict
      if !haskey(gkm.w, e) || gkm.w[e] != w_e
        printDiagnostics && println("Edge $e weight inconsistency in w dict")
        return false
      end
    end
  end

  if length(gkm.labels) != nv
    printDiagnostics && println("Not the right number of labels")
    return false
  elseif (val > 1 && !is2_indep(gkm))
    printDiagnostics && println("GKM graph is not 2-independent.")
    return false
  end

  if length(unique(gkm.labels)) != length(gkm.labels)
    printDiagnostics && println("Labels are not unique.")
    return false
  end

  # The output below is more annoying than useful for dealing with non-3-independent GKM graphs which are perfectly fine.

  # if (val > 2 && !is3_indep(gkm))
  #   printDiagnostics && println("GKM graph is valid but not 3-independent, so connections may not be unique.")
  # end

  if length(gens(gkm.equivariantCohomology.coeffRing)) != rank_torus(gkm)
    printDiagnostics && println("Coefficient ring of equivariant cohomology has wrong rank.")
    return false
  end

  if length(gens(gkm.equivariantCohomology.cohomRing)) != nv
    printDiagnostics && println("Equivariant cohomology ring should have rank = number of vertices.")
    return false
  end

  return true
end

function edgeFromLabels(G::AbstractGKM_graph, s::String, d::String)::Edge
  @req s in G.labels "source $s is not a vertex in $G."
  @req d in G.labels "destination $d is not a vertex in $G."

  sd = indexin([s, d], G.labels)
  return Edge(sd[1], sd[2])
end

@doc raw"""
    enlarge_torus(G::AbstractGKM_graph, r::Int64) -> AbstractGKM_graph

Return a copy of the GKM graph $G$ where the acting torus has been enlarged
by $r$ dimensions which act trivially.
This function is sometimes used in the construction of vector bundles, when one would like
the vector bundle and the base space to have the same weight lattice.

# Example
```jldoctest enlarge_torus_test
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> enlarge_torus(P2, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0, 0, 0)
3 -> 1 => (-1, 0, 1, 0, 0)
3 -> 2 => (0, -1, 1, 0, 0)
```
"""
function enlarge_torus(G::AbstractGKM_graph, r::Int64)::AbstractGKM_graph
  @req r>=0 "r must be non-negative"
  r == 0 && return G

  r1 = rank_torus(G)
  weightType = parent(_get_weight_type(G))
  M2 = free_module(weightType, r1 + r)
  g2 = gens(M2)
  f_inc = ModuleHomomorphism(G.M, M2, [g2[i] for i in 1:r1]);
  return substitute_torus(G, f_inc)
end

@doc raw"""
    substitute_torus(G::AbstractGKM_graph{R}, f::AbstractAlgebra.Generic.ModuleHomomorphism{R}) where R <: GKM_weight_type

Return a copy of the GKM graph `G` where the weights of all flags and edges are substituted
according to the module homomorphism `f`.

If `G` has a natural connection (i.e. if it is 3-independent or a connection has been set using `set_connection!`),
then it is copied to the result.

# Example

```jldoctest substitute_torus
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> M = free_module(ZZ, 2);

julia> g1, g2 = gens(M);

julia> f = ModuleHomomorphism(G.M, M, [zero(M), g1, g2]);

julia> G2 = substitute_torus(G, f)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0)
3 -> 1 => (0, 1)
3 -> 2 => (-1, 1)
```
"""
function substitute_torus(G::AbstractGKM_graph{R}, f::AbstractAlgebra.Generic.ModuleHomomorphism{R}) where R <: GKM_weight_type
  M_new = codomain(f)
  M_old = domain(f)
  @req G.M == M_old "Domain of f must be G.M"

  nv = n_vertices(G.g)
  weights_at_vertex = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}}(undef, nv)
  w = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}}()

  for i in 1:nv
    weights_at_vertex[i] = [f(w) for w in G.weights_at_vertex[i]]
  end

  for e in keys(G.w)
    w[e] = f(G.w[e])
  end

  # create deepcopy so that later modifications to G don't affect the result.
  GW_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()
  res = AbstractGKM_graph(deepcopy(G.g), deepcopy(G.labels), M_new, weights_at_vertex, deepcopy(G.edge_to_flag_index), deepcopy(G.flag_to_edge), w, nothing, nothing, deepcopy(G.connection), GW_structure_consts, false)
  if !isnothing(res.connection)
    res.connection.gkm = res
  end
  res.equivariantCohomology = _equivariant_cohomology_ring(res)

  return res
end