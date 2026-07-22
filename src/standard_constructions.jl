export flag_variety, grassmannian, gkm_graph_of_toric, projective_space, schubert_class, schubert_classes

function _type_a_companion_root(
  alpha::NTuple{2,Int},
  beta::NTuple{2,Int},
)::Union{Nothing,NTuple{2,Int}}
  i, j = alpha
  k, l = beta

  alpha == beta && return nothing

  if i == k
    return (j, l)
  elseif j == l
    return (k, i)
  elseif j == k
    return (i, l)
  elseif l == i
    return (k, j)
  end

  return nothing
end

function _type_a_geometric_connection_integer(
  alpha::NTuple{2,Int},
  beta::NTuple{2,Int},
  block_of_position::Vector{Int},
)::ZZRingElem
  alpha == beta && return ZZ(2)

  gamma = _type_a_companion_root(alpha, beta)
  gamma === nothing && return ZZ(0)

  p, q = gamma
  is_positive_omitted =
    p < q && block_of_position[p] != block_of_position[q]

  return is_positive_omitted ? ZZ(0) : ZZ(1)
end

function _type_a_combinatorial_connection_integer(
  alpha::NTuple{2,Int},
  beta::NTuple{2,Int},
)::ZZRingElem
  i, j = alpha
  k, l = beta

  return ZZ(
    Int(i == k) + Int(j == l) -
    Int(i == l) - Int(j == k),
  )
end

@doc raw"""
    flag_variety(::Type{GKM_graph}, s::Vector{Int64}; connection::Symbol = :geometric) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the variety of flags of ``\mathbb{C}^n``. The dimensions of quotients are expressed by the array `s`. The labels represent the vectors generating the flags. For example, if ``s=[1,2,1]``, the string ``213`` corresponds to the flag:

``0\subset \langle e_2 \rangle \subset \langle e_2, e_1, e_3 \rangle \subset \langle e_2, e_1, e_3, e_4 \rangle=\mathbb{C}^4.``

!!! note
    This function is faster than `generalized_gkm_flag(root_system(:A, n-1), S)`, but the results are isomorphic.

# Choice of compatible connection

The optional argument `connection` selects the compatible connection stored on the graph.
The options are the same as in the more general [`generalized_gkm_flag`](@ref).

- `:geometric` (default): uses the Birkhoff-Grothendieck splitting of the tangent bundle along
  each torus invariant ``\mathbb{P}^1``, following [McKay_Benjamin_2006; Lemma 16](@cite).
- `:combinatorial`: uses the root-label-preserving connection described in
  [Guillemin_Holm_Zara_2006; Section 2.2.7](@cite).

# Examples
```jldoctest
julia> flag_variety(GKM_graph, [1,3])
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> flag_variety(GKM_graph, [2,1])
GKM graph with 3 nodes, valency 2 and axial function:
13 -> 12 => (0, -1, 1)
23 -> 12 => (-1, 0, 1)
23 -> 13 => (-1, 1, 0)

```

Let us also see the difference between the geometric and the combinatorial connection on the full flag variety of $\mathbb{C}^3$.

```jldoctest
julia> F3_with_geometric_con = flag_variety(GKM_graph, [1,1,1]; connection=:geometric)
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

julia> F3_with_combinatorial_con = flag_variety(GKM_graph, [1,1,1]; connection=:combinatorial)
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

julia> C_combinatorial = get_connection(F3_with_combinatorial_con)
GKM connection for GKM graph with 6 nodes and valency 3:
Connection:
Edge(4, 2) => [2, 1, 3]
Edge(3, 5) => [1, 3, 2]
Edge(3, 1) => [2, 1, 3]
Edge(4, 3) => [3, 2, 1]
Edge(5, 2) => [3, 2, 1]
Edge(1, 3) => [2, 1, 3]
Edge(6, 5) => [2, 1, 3]
Edge(1, 2) => [1, 3, 2]
Edge(5, 3) => [1, 3, 2]
Edge(2, 4) => [2, 1, 3]
Edge(6, 1) => [3, 2, 1]
Edge(2, 1) => [1, 3, 2]
Edge(4, 6) => [1, 3, 2]
Edge(3, 4) => [3, 2, 1]
Edge(1, 6) => [3, 2, 1]
Edge(5, 6) => [2, 1, 3]
Edge(6, 4) => [1, 3, 2]
Edge(2, 5) => [3, 2, 1]
a_i's:
Edge(4, 2) => ZZRingElem[2, 1, 1]
Edge(3, 5) => ZZRingElem[1, 1, 2]
Edge(3, 1) => ZZRingElem[2, -1, 1]
Edge(4, 3) => ZZRingElem[1, 2, -1]
Edge(5, 2) => ZZRingElem[2, 1, -1]
Edge(1, 3) => ZZRingElem[-1, 2, 1]
Edge(6, 5) => ZZRingElem[1, -1, 2]
Edge(1, 2) => ZZRingElem[2, -1, 1]
Edge(5, 3) => ZZRingElem[1, 2, 1]
Edge(2, 4) => ZZRingElem[1, 2, 1]
Edge(6, 1) => ZZRingElem[2, 1, 1]
Edge(2, 1) => ZZRingElem[2, 1, -1]
Edge(4, 6) => ZZRingElem[1, -1, 2]
Edge(3, 4) => ZZRingElem[-1, 2, 1]
Edge(1, 6) => ZZRingElem[1, 1, 2]
Edge(5, 6) => ZZRingElem[-1, 1, 2]
Edge(6, 4) => ZZRingElem[1, 2, -1]
Edge(2, 5) => ZZRingElem[-1, 1, 2]

julia> C_geometric = get_connection(F3_with_geometric_con)
GKM connection for GKM graph with 6 nodes and valency 3:
Connection:
Edge(4, 2) => [2, 1, 3]
Edge(3, 5) => [1, 3, 2]
Edge(3, 1) => [2, 3, 1]
Edge(4, 3) => [1, 2, 3]
Edge(5, 2) => [3, 1, 2]
Edge(1, 3) => [3, 1, 2]
Edge(6, 5) => [1, 2, 3]
Edge(1, 2) => [1, 2, 3]
Edge(5, 3) => [1, 3, 2]
Edge(2, 4) => [2, 1, 3]
Edge(6, 1) => [3, 2, 1]
Edge(2, 1) => [1, 2, 3]
Edge(4, 6) => [3, 1, 2]
Edge(3, 4) => [1, 2, 3]
Edge(1, 6) => [3, 2, 1]
Edge(5, 6) => [1, 2, 3]
Edge(6, 4) => [2, 3, 1]
Edge(2, 5) => [2, 3, 1]
a_i's:
Edge(4, 2) => ZZRingElem[2, 1, 1]
Edge(3, 5) => ZZRingElem[1, 1, 2]
Edge(3, 1) => ZZRingElem[2, 0, 0]
Edge(4, 3) => ZZRingElem[0, 2, 0]
Edge(5, 2) => ZZRingElem[2, 0, 0]
Edge(1, 3) => ZZRingElem[0, 2, 0]
Edge(6, 5) => ZZRingElem[0, 0, 2]
Edge(1, 2) => ZZRingElem[2, 0, 0]
Edge(5, 3) => ZZRingElem[1, 2, 1]
Edge(2, 4) => ZZRingElem[1, 2, 1]
Edge(6, 1) => ZZRingElem[2, 1, 1]
Edge(2, 1) => ZZRingElem[2, 0, 0]
Edge(4, 6) => ZZRingElem[0, 0, 2]
Edge(3, 4) => ZZRingElem[0, 2, 0]
Edge(1, 6) => ZZRingElem[1, 1, 2]
Edge(5, 6) => ZZRingElem[0, 0, 2]
Edge(6, 4) => ZZRingElem[0, 2, 0]
Edge(2, 5) => ZZRingElem[0, 0, 2]
```
"""
function flag_variety(
  ::Type{GKM_graph},
  s::Vector{Int64};
  connection::Symbol=:geometric,
)
  _validate_homogeneous_connection_option(connection)

  @req !isempty(s) "the vector of dimensions is empty"
  @req all(i -> s[i] > 0, eachindex(s)) "all dimensions must be positive"

  K::Vector{Int64} = [sum(s[1:i]) for i in 0:length(s)]
  n = K[end]

  block_of_position = Vector{Int}(undef, n)
  for b in eachindex(s)
    block_of_position[(K[b] + 1):K[b + 1]] .= b
  end

  omitted_roots = NTuple{2,Int}[
    (i, j) for i in 1:(n - 1) for j in (i + 1):n if
    block_of_position[i] != block_of_position[j]
  ]
  omitted_root_set = Set(omitted_roots)

  a_by_root_pair = Dict{
    Tuple{NTuple{2,Int},NTuple{2,Int}},
    ZZRingElem,
  }()
  for alpha in omitted_roots
    for beta in omitted_roots
      a_by_root_pair[(alpha, beta)] = if connection === :geometric
        _type_a_geometric_connection_integer(alpha, beta, block_of_position)
      else
        @assert connection === :combinatorial
        _type_a_combinatorial_connection_integer(alpha, beta)
      end
    end
  end

  d::Dict{Int64, NTuple{K[end], Int64}} = Dict{Int64, NTuple{K[end], Int64}}()
  index::Int64 = 1
  
  for c in Combinatorics.permutations(1:K[end])
    if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:length(s))
      d[index] = (c...,)
      index += 1
    end
  end
      
  nv::Int64 = length(keys(d))
  g = Graph{Undirected}(nv)
  position_of_value = [invperm(collect(d[v])) for v in 1:nv]

  M = free_module(ZZ, K[end])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  local_root = Dict{Edge,NTuple{2,Int}}()

  for v in 1:nv
    for w in (v+1):nv
      dif = _T_link(d[v], d[w], K)
        
      isempty(dif) && continue

      add_edge!(g, v, w)
      W[Edge(w, v)] = gens(M)[dif[2]] - gens(M)[dif[1]]

      x, y = dif
      root_at_v = minmax(position_of_value[v][x], position_of_value[v][y])
      root_at_w = minmax(position_of_value[w][x], position_of_value[w][y])

      @assert root_at_v in omitted_root_set
      @assert root_at_w in omitted_root_set

      local_root[Edge(v, w)] = root_at_v
      local_root[Edge(w, v)] = root_at_w
    end
  end

  @assert all(
    haskey(local_root, e) && haskey(local_root, reverse(e)) for e in edges(g)
  ) "Every directed type A flag edge must have a local positive-root label"

  a = Dict{Tuple{Edge,Edge},ZZRingElem}()
  for v in vertices(g)
    outgoing = [Edge(v, u) for u in all_neighbors(g, v)]
    for e in outgoing
      alpha = local_root[e]
      for e_prime in outgoing
        beta = local_root[e_prime]
        a[(e, e_prime)] = a_by_root_pair[(alpha, beta)]
      end
    end
  end

  labels = [prod(i -> d[n][i]>9 ? "|$(d[n][i])|" : "$(d[n][i])", 1:K[end-1]) for n in 1:nv] 

  G = gkm_graph(g, labels, M, W)


  con = build_GKM_connection(G, a)
  @req isvalid(con; printDiagnostics=false) "Invalid $(connection) connection for type A flag variety"
  set_connection!(G, con)

  return G

end

  
# this function variefies if two flags can be connected by a T-invariant curve. If this is not possible, an empty array is returned. Otherwise, it returns the two indices that need to be swapped
function _T_link(a::Tuple{Vararg{Int64}}, b::Tuple{Vararg{Int64}}, K::Vector{Int64})
    
  ans::Vector{Int64} = [0, 0]
  
  for i in 1:length(K)
    a1 = a[1:K[i]]
    b1 = b[1:K[i]]
    a1_b1 = setdiff(a1, b1)

    if isempty(a1_b1)
      continue
    end

    if length(a1_b1) > 1
      return Int64[]
    end

    b1_a1 = setdiff(b1, a1)

    if ans[1] == 0
      ans = [a1_b1[1], b1_a1[1]]
      continue
    end

    if ans != [a1_b1[1], b1_a1[1]]
      return Int64[]
    end
end

  return sort(ans)
end

@doc raw"""
    grassmannian(::Type{GKM_graph}, k::Int, n::Int; connection::Symbol = :geometric) -> AbstractGKM_graph{ZZRingElem}

Construct the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

The optional argument `connection` accepts `:geometric` (default) and
`:combinatorial` with the same meaning as in [`flag_variety`](@ref).

# Examples
```jldoctest
julia> grassmannian(GKM_graph, 2, 4)
GKM graph with 6 nodes, valency 4 and axial function:
13 -> 12 => (0, -1, 1, 0)
14 -> 12 => (0, -1, 0, 1)
14 -> 13 => (0, 0, -1, 1)
23 -> 12 => (-1, 0, 1, 0)
23 -> 13 => (-1, 1, 0, 0)
24 -> 12 => (-1, 0, 0, 1)
24 -> 14 => (-1, 1, 0, 0)
24 -> 23 => (0, 0, -1, 1)
34 -> 13 => (-1, 0, 0, 1)
34 -> 14 => (-1, 0, 1, 0)
34 -> 23 => (0, -1, 0, 1)
34 -> 24 => (0, -1, 1, 0)

```
"""
function grassmannian(
  ::Type{GKM_graph},
  k::Int,
  n::Int;
  connection::Symbol=:geometric,
)
  _validate_homogeneous_connection_option(connection)
  @req (k >= 0 && n >= k) "Dimension must be non-negative"
  
  return flag_variety(GKM_graph, [k, n-k]; connection=connection)
end

@doc raw"""
    projective_space(::Type{gkm_graph}, d::Int)

Construct the projective space of dimension `d`.

# Examples
```jldoctest
julia> projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
```
"""
function projective_space(::Type{GKM_graph}, d::Int)
  @req d >= 0 "Dimension must be non-negative"
  
  return grassmannian(GKM_graph, 1, d+1)
end


@doc raw"""
    gkm_graph_of_toric(v::NormalToricVariety; small_torus::Bool=false) -> AbstractGKM_graph{ZZRingElem}

Construct the GKM graph of the smooth toric variety `v`.

If the variety is projective, all flags will be connected to edges.
If the variety is non-projective, codimension-1 faces that are not shared by two maximal cones
will correspond to standalone flags.

# Dimension of the torus
If the optional argument `small_torus` is `false` (default value) then the torus rank of the
result is the number of rays of `v`.
If `small_torus` is `true` then the torus rank of the result is the dimension of `v`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> gkm_graph_of_toric(P2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1)
3 -> 1 => (0, 1, -1)
3 -> 2 => (-1, 1, 0)

julia> gkm_graph_of_toric(P2; small_torus=true)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (1, 0)
3 -> 1 => (0, 1)
3 -> 2 => (-1, 1)

julia> F = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> gkm_graph_of_toric(F)
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0, -1, 0)
3 -> 2 => (3, 1, 0, -1)
4 -> 1 => (0, 1, 3, -1)
4 -> 3 => (-1, 0, 1, 0)

julia> gkm_graph_of_toric(F; small_torus=true)
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (1, 0)
3 -> 2 => (3, 1)
4 -> 1 => (0, 1)
4 -> 3 => (-1, 0)

julia> gkm_graph_of_toric(affine_space(NormalToricVariety, 3))
GKM graph with 1 nodes, valency 3 and axial function:
Standalone flags:
1.1 => (-1, 0, 0)
1.2 => (0, -1, 0)
1.3 => (0, 0, -1)
```
"""
function gkm_graph_of_toric(v::NormalToricVariety; small_torus::Bool=false)

  @req is_smooth(v) "toric variety must be smooth"

  len = length(maximal_cones(v))
  g = Graph{Undirected}(len)
  M = free_module(ZZ, small_torus ? dim(v) : n_rays(v))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()

  # Build the graph with edges for shared codimension-1 faces
  for sigma1 in 1:(len-1)
    for sigma2 in (sigma1+1):len
      sigma1_cone = maximal_cones(v)[sigma1]
      sigma2_cone = maximal_cones(v)[sigma2]

      # Check if these cones share a codimension-1 face
      count(x -> x in rays(sigma1_cone), rays(sigma2_cone)) != (dim(v) - 1) && continue

      # Find the ray in sigma1 that's not in sigma2
      ray1 = findfirst(r -> !(r in rays(sigma2_cone)), rays(sigma1_cone))

      add_edge!(g, sigma1, sigma2)
      W[Edge(sigma2, sigma1)] = _omega(v, sigma1, ray1, M; small_torus)
    end
  end

  # Create the base GKM graph
  G = gkm_graph(g, ["$i" for i in 1:len], M, W; check=false) # don't check as valency is wrong until flags are added.

  # Now add standalone flags for codimension-1 faces that don't belong to another maximal cone
  for sigma_idx in 1:len
    sigma = maximal_cones(v)[sigma_idx]

    # For each ray in this cone
    for (ray_idx, ray) in enumerate(rays(sigma))
      # Check if this ray corresponds to a codimension-1 face shared with another cone
      # by checking if removing this ray gives a face contained in another maximal cone
      codim1_face = [r for r in rays(sigma) if r != ray]

      # Check if this codim-1 face is contained in any other maximal cone
      is_shared = false
      for other_idx in 1:len
        other_idx == sigma_idx && continue

        other_sigma = maximal_cones(v)[other_idx]
        if all(r -> r in rays(other_sigma), codim1_face)
          is_shared = true
          break
        end
      end

      # If not shared, this is a standalone flag
      if !is_shared
        weight = _omega(v, sigma_idx, ray_idx, M; small_torus)
        add_standalone_flag!(G, sigma_idx, -weight)
      end
    end
  end

  @req isvalid(G) "gkm_graph_of_toric produced invalid result."

  return G
end

function _omega(v::NormalToricVariety, n_SIGMA::Int64, ray_idx::Int64, M::AbstractAlgebra.Generic.FreeModule{ZZRingElem}; small_torus::Bool=false)::AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}

  SIGMA = maximal_cones(v)[n_SIGMA]
  ray = rays(SIGMA)[ray_idx]
  scalars = gens(M)
  ud = QQFieldElem[]

  # Find a polarizing ray that pairs non-trivially with our ray
  for pol_ray in rays(polarize(SIGMA))
    if dot(ray, pol_ray) != 0
      ud = lcm(denominator.(pol_ray)) * pol_ray
      break
    end
  end

  ans = zero(M)

  if small_torus
    for k in 1:dim(v)
      ans += Int64(ud[k]) * scalars[k]
    end
  else
    for (k, vi) in enumerate(rays(v))
      ans += Int64(dot(vi, ud)) * scalars[k]
    end
  end

  return ans
end
