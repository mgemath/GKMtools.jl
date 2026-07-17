@doc raw"""
    flag_variety(::Type{GKMGraph}, s::AbstractVector{<:Integer}) -> GKMGraph

Construct the GKM graph of the partial flag variety whose successive quotient
dimensions are given by `s`.

The implementation generates only the required block permutations and builds
edges by directly exchanging entries between distinct blocks.
The labels represent the vectors generating the flags. For example, if ``s=[1,2,1]``, the string ``213`` corresponds to the flag:

``0\subset \langle e_2 \rangle \subset \langle e_2, e_1, e_3 \rangle \subset \langle e_2, e_1, e_3, e_4 \rangle=\mathbb{C}^4.``

!!! note
    This function is faster than `generalized_gkm_flag(root_system(:A, n-1), S)`, but the results are isomorphic.

# Examples
```jldoctest
julia> flag_variety(GKMGraph, [1,3])
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)
Homogeneous connection for GKM graph with 4 nodes and valency 3

julia> flag_variety(GKMGraph, [2,1])
GKM graph with 3 nodes, valency 2 and axial function:
13 -> 12 => (0, -1, 1)
23 -> 12 => (-1, 0, 1)
23 -> 13 => (-1, 1, 0)
Homogeneous connection for GKM graph with 3 nodes and valency 2
```
"""
function flag_variety(
  ::Type{GKMGraph},
  s::AbstractVector{<:Integer},
)
  @req !isempty(s) "the vector of dimensions is empty"
  @req all(>(0), s) "all dimensions must be positive"

  return _flag_variety_graph(Int.(s))
end


function _flag_variety_graph(s::Vector{Int})
  n = sum(s)
  cuts = cumsum(vcat(0, s))

  ###########################################################################
  # Vertices
  ###########################################################################

  # A multiset permutation records the block containing each integer
  # 1, ..., n. This generates exactly n! / prod(s_i!) representatives,
  # instead of generating all n! permutations and filtering afterward.
  block_data = reduce(
    vcat,
    (fill(block, s[block]) for block in eachindex(s)),
  )

  block_words = Combinatorics.multiset_permutations(block_data, n)



  representatives = Tuple[
    Tuple(
      value
      for block in eachindex(s)
      for value in 1:n
      if word[value] == block
    )
    for word in block_words
  ]

  vertex_index = Dict(
    representative => index
    for (index, representative) in enumerate(representatives)
  )
  flag_type = NTuple{length(representatives[1]), Int}
  ###########################################################################
  # Graph and roots
  ###########################################################################

  g = Graph{Undirected}(length(representatives))

  # Store the two coordinate indices associated with each unoriented edge.
  roots = Dict{Tuple{Int,Int},Tuple{Int,Int}}()

  if length(s) > 1
    for v in eachindex(representatives)
      representative = collect(representatives[v])

      for first_block in 1:(length(s) - 1)
        first_range =
          (cuts[first_block] + 1):cuts[first_block + 1]

        for second_block in (first_block + 1):length(s)
          second_range =
            (cuts[second_block] + 1):cuts[second_block + 1]

          for p in first_range
            for q in second_range
              swapped = copy(representative)
              swapped[p], swapped[q] = swapped[q], swapped[p]

              # Restore the canonical representative: entries within every
              # block are sorted.
              sort!(@view swapped[first_range])
              sort!(@view swapped[second_range])

              w = vertex_index[Tuple(swapped)]

              # Each edge is encountered from both endpoints.
              v < w || continue

              add_edge!(g, v, w)
              roots[(v, w)] =
                minmax(representative[p], representative[q])
            end
          end
        end
      end
    end
  end

  ###########################################################################
  # Flags and axial weights
  ###########################################################################

  M = free_module(ZZ, n)
  basis = gens(M)

  flags = [
    FlagWeight{ZZRingElem}[]
    for _ in vertices(g)
  ]

  edge_flags = Dict{Edge,Tuple{Int,Int}}()

  for e in edges(g)
    v = src(e)
    w = dst(e)

    low, high = roots[minmax(v, w)]
    canonical_weight = basis[high] - basis[low]

    # OSCAR's canonical orientation for these edges points from the larger
    # vertex index to the smaller one.
    weight_at_v =
      v > w ? canonical_weight : -canonical_weight

    push!(
      flags[v],
      FlagWeight{ZZRingElem}(weight_at_v),
    )

    push!(
      flags[w],
      FlagWeight{ZZRingElem}(-weight_at_v),
    )

    edge_flags[e] = (
      length(flags[v]),
      length(flags[w]),
    )
  end

  ###########################################################################
  # Typed vertices and new GKM graph
  ###########################################################################

  visible_length = cuts[end - 1]

  labels = FlagVertex[
    FlagVertex{flag_type}(_flag_label(representative, visible_length), representative)
    for representative in representatives
  ]

  core = GKMCombinatorialData{
    ZZRingElem,
    FlagVertex{flag_type},
    FlagWeight{ZZRingElem},
  }(
    g,
    M,
    labels,
    flags,
    edge_flags,
  )

  connection = connection = _homogeneous_connection(core, roots)
  cohomology = create_cohomology(n, length(representatives))

  return GKMGraph{
    ZZRingElem,
    FlagVertex{flag_type},
    FlagWeight{ZZRingElem},
  }(
    core,
    connection,
    cohomology,
    nothing,
    nothing,
  )
end


function _flag_label(
  representative::Tuple,
  visible_length::Int,
)
  visible_length == 0 && return "1"

  return join(
    value > 9 ? "|$value|" : string(value)
    for value in representative[1:visible_length]
  )
end


@doc raw"""
    grassmannian(::Type{GKMGraph}, k::Integer, n::Integer) -> GKMGraph

Construct the GKM graph of the Grassmannian of `k`-planes in
``\mathbb{C}^n``.
# Examples
```jldoctest
julia> grassmannian(GKMGraph, 2, 4)
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
Homogeneous connection for GKM graph with 6 nodes and valency 4
```
"""
function grassmannian(
  ::Type{GKMGraph},
  k::Integer,
  n::Integer,
)
  @req 0 < k < n "require 0 < k < n"

  return _flag_variety_graph([
    Int(k),
    Int(n - k),
  ])
end


@doc raw"""
    projective_space(::Type{GKMGraph}, d::Integer) -> GKMGraph

Construct complex projective space of dimension `d`.
```jldoctest
julia> projective_space(GKMGraph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Homogeneous connection for GKM graph with 3 nodes and valency 2
"""
function projective_space(
  ::Type{GKMGraph},
  d::Integer,
)
  @req d > 0 "dimension must be positive"

  return grassmannian(GKMGraph, 1, d + 1)
end

function _homogeneous_connection(
  core::GKMCombinatorialData{
    ZZRingElem,
    V,
    FlagWeight{ZZRingElem},
  },
  roots::Dict{Tuple{Int,Int},Tuple{Int,Int}},
) where {V<:FlagVertex}
  transport = Dict{Edge,Vector{Int}}()
  coefficients = Dict{Edge,Vector{ZZRingElem}}()

  # Record the root associated with each flag at each vertex.
  flag_roots = [
    Vector{Tuple{Int,Int}}(undef, length(flags(core, v)))
    for v in vertices(core)
  ]

  for e in edges(graph(core))
    v = src(e)
    w = dst(e)

    root = roots[minmax(v, w)]
    source_flag, target_flag = core.edge_flags[e]

    flag_roots[v][source_flag] = root
    flag_roots[w][target_flag] = root
  end

  # Root-to-flag lookup at each vertex makes construction of the transport
  # linear in the number of flags.
  root_to_flag = [
    Dict(root => i for (i, root) in enumerate(flag_roots[v]))
    for v in vertices(core)
  ]

  for unoriented_edge in edges(graph(core))
    for e in (unoriented_edge, reverse(unoriented_edge))
      v = src(e)
      w = dst(e)

      low, high = roots[minmax(v, w)]
      edge_weight = weight(core, e)

      image = Vector{Int}(undef, length(flag_roots[v]))
      edge_coefficients =
        Vector{ZZRingElem}(undef, length(flag_roots[v]))

      for i in eachindex(flag_roots[v])
        beta_low, beta_high = flag_roots[v][i]

        # Transport is induced by the transposition corresponding to the
        # invariant curve.
        transported_root = minmax(
          _apply_transposition(beta_low, low, high),
          _apply_transposition(beta_high, low, high),
        )

        j = root_to_flag[w][transported_root]
        image[i] = j

        difference =
          weight(core, v, i) -
          weight(core, w, j)

        edge_coefficients[i] =
          _homogeneous_connection_coefficient(
            difference,
            edge_weight,
          )
      end

      transport[e] = image
      coefficients[e] = edge_coefficients
    end
  end

  return Connection{ZZRingElem}(
    transport,
    coefficients,
    "Homogeneous",
  )
end


function _apply_transposition(
  value::Int,
  first::Int,
  second::Int,
)
  value == first && return second
  value == second && return first
  return value
end


function _homogeneous_connection_coefficient(
  difference,
  edge_weight,
)
  coordinate = findfirst(
    i -> !iszero(edge_weight[i]),
    1:rank(parent(edge_weight)),
  )

  @assert !isnothing(coordinate) "edge weight must be nonzero"

  coefficient = divexact(
    difference[coordinate],
    edge_weight[coordinate],
  )

  @assert difference == coefficient * edge_weight """
  incompatible homogeneous connection
  """

  return coefficient
end