function _bundle_weight_class(
  E::AbstractGKMVectorBundle,
  v::Int,
  i::Int,
  t::Vector{T},
) where {T}
  G = baseof(E)

  @req rank(E.M) == rank_torus(G) """
  Chern classes are currently supported when bundle and base have the same torus rank
  """

  res = zero(t[1])
  w = fiber_weight(E, v, i)

  for j in 1:rank_torus(G)
    res += w[j] * t[j]
  end

  return res
end

function _elementary_symmetric_class(values, k::Int)
  @req k >= 0 """
  Chern class is only defined for non-negative index
  """

  k == 0 && return one(first(values))
  k > length(values) && return zero(first(values))

  res = zero(first(values))

  for combination in Combinatorics.combinations(
    eachindex(values),
    k,
  )
    res += prod(
      (values[i] for i in combination);
      init=one(first(values)),
    )
  end

  return res
end

@doc raw"""
    chern_class(E::AbstractGKMVectorBundle, k::Integer)
    chern_class(G::AbstractGKMGraph, k::Integer)

Return the `k`-th equivariant Chern class of a smooth or orbifold GKM vector
bundle `E`. For a graph `G`, return the Chern class of its tangent bundle.

The result is a localized [`GKMClass`](@ref), whose restriction at each vertex
is the `k`-th elementary symmetric polynomial in the fiber weights. The zeroth
Chern class is the unit, and classes above the bundle rank are zero.

Throw an error when `k` is negative or the torus ranks of the bundle and its
base differ.

# Examples
In the first example, we calculate the Chern classes of $\mathbb{P}^2$ and its Chern numbers $c_1^2$ and $c_2$.
```jldoctest chern_class_P2
julia> G = projective_space(GKMGraph, 2);

julia> chern_class(G, 0)
GKM class with restrictions: 
[1, 1, 1]

julia> chern_class(G, 1)
GKM class with restrictions:
[2*t1 - t2 - t3, -t1 + 2*t2 - t3, -t1 - t2 + 2*t3]

julia> chern_class(G, 2)
GKM class with restrictions:
[t1^2 - t1*t2 - t1*t3 + t2*t3, -t1*t2 + t1*t3 + t2^2 - t2*t3, t1*t2 - t1*t3 - t2*t3 + t3^2]

julia> chern_class(G, 3)
GKM class with restrictions:
[0, 0, 0]

julia> integrate(chern_class(G, 1)*chern_class(G, 1))
9

julia> integrate(chern_class(G, 2))
3
``` 
In the second example, we calculate the Chern number $c_1c_2$ for the twisted flag variety (see [`gkm_3d_twisted_flag`](@ref GKMtools.gkm_3d_twisted_flag)).
This number coincides with the sum of the Chern numbers of the edges by [GS14; Proposition 4.6](@cite).
```jldoctest chern_class_twisted_flag
julia> G = gkm_3d_twisted_flag();

julia> integrate(chern_class(G, 1) * chern_class(G, 2))
24

julia> sum([chern_number(e, G) for e in edges(G)])
24
```


Let us compute the chern classes of the equivariant vector bundle
$\mathcal{O}(-2)\oplus\mathcal{O}(3)\rightarrow\mathbb{P}^2$ where each summand has its own
equivariant parameter.

```jldoctest vec_bdle_chern_test
julia> O1 = vector_bundle_O(2, [1])

julia> V = O1^(-2)+O1^(3);

julia> c1 = chern_class(V, 1)
GKM class with restrictions: 
[t1, t2, t3]

julia> c2 = chern_class(V, 2)
GKM class with restrictions: 
[-6*t1^2, -6*t2^2, -6*t3^2]

julia> integrate(c2)
-6
```
"""
function Oscar.chern_class(
  E::AbstractGKMVectorBundle,
  k::Integer,
)
  @req k >= 0 """
  Chern class is only defined for non-negative index
  """

  G = baseof(E)
  basis = gens_cohomRing(G)

  k == 0 && return one(first(basis))
  k > rank(E) && return zero(first(basis))

  parameters = gens_coeffRing(G)
  result = zero(first(basis))

  for v in 1:num_vertices(G)
    local_weights = [
      _bundle_weight_class(E, v, i, parameters)
      for i in 1:rank(E)
    ]

    result +=
      _elementary_symmetric_class(
        local_weights,
        Int(k),
      ) * basis[v]
  end

  return result
end

@doc raw"""
    first_chern_class(E::AbstractGKMVectorBundle)
    first_chern_class(G::AbstractGKMGraph)

Return the first equivariant Chern class of `E`, or of the tangent bundle of
`G`. This is equivalent to `chern_class(E, 1)` or `chern_class(G, 1)`.


# Examples
```jldoctest first_chern_class
julia> P2 = projective_space(GKM_graph, 2);

julia> first_chern_class(P2)
GKM class with restrictions: 
[2*t1 - t2 - t3, -t1 + 2*t2 - t3, -t1 - t2 + 2*t3]

julia> H3 = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 3))
GKM graph with 4 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (-3, -1, 0, 1)
4 -> 1 => (0, -1, -3, 1)
4 -> 3 => (1, 0, -1, 0)
Algorithmic connection for GKM graph with 4 nodes and valency 2

julia> first_chern_class(H3)
GKM class with restrictions: 
[t1 + t2 + 2*t3 - t4, 2*t1 + t2 + t3 - t4, -4*t1 - t2 + t3 + t4, t1 - t2 - 4*t3 + t4]
```
"""
first_chern_class(
  E::AbstractGKMVectorBundle,
) = chern_class(E, 1)

"""
    chern_classes(E::AbstractGKMVectorBundle)
    chern_classes(G::AbstractGKMGraph)

Return all equivariant Chern classes from degree zero through the rank. For a
graph, compute the classes of its tangent bundle. The first element is the
unit class.

# Examples
```jldoctest
julia> X = hirzebruch_surface(NormalToricVariety, 2);

julia> E = [toric_line_bundle(X, [1, 0]), toric_line_bundle(X, [0, 1])];

julia> V = gkm_vector_bundle_of_toric(E);

julia> (c0, c1, c2) = chern_classes(V);

julia> c0
GKM class with restrictions: 
[1, 1, 1, 1]

julia> c1
GKM class with restrictions: 
[-t3 - t4 - t5 + t6, -t1 - t3 - t4 + t6, t1 + t2 - t3 - t4, t2 - t3 - t4 + t5]

julia> c2
GKM class with restrictions: 
[t3*t4 + 2*t3*t5 - t3*t6 - t4*t5 - 2*t5^2 + t5*t6, -2*t1^2 + 2*t1*t3 - t1*t4 + t1*t6 + t3*t4 - t3*t6, t1*t2 - t1*t4 - t2*t3 + t3*t4, -t2*t3 + t2*t5 + t3*t4 - t4*t5]
```
"""
function Oscar.chern_classes(
  E::AbstractGKMVectorBundle,
)
  return [
    chern_class(E, k)
    for k in 0:rank(E)
  ]
end

@doc raw"""
    total_chern_class(E::AbstractGKMVectorBundle)
    total_chern_class(G::AbstractGKMGraph)

Return the total equivariant Chern class, defined as the sum of all Chern
classes of `E`. For a graph, compute the total Chern class of its tangent
bundle.

# Examples
```jldoctest
julia> X = hirzebruch_surface(NormalToricVariety, 2);

julia> E = [toric_line_bundle(X, [1, 0]), toric_line_bundle(X, [0, 1])];

julia> V = gkm_vector_bundle_of_toric(E);

julia> total_chern_class(V)
GKM class with restrictions: 
[t3*t4 + 2*t3*t5 - t3*t6 - t3 - t4*t5 - t4 - 2*t5^2 + t5*t6 - t5 + t6 + 1, -2*t1^2 + 2*t1*t3 - t1*t4 + t1*t6 - t1 + t3*t4 - t3*t6 - t3 - t4 + t6 + 1, t1*t2 - t1*t4 + t1 - t2*t3 + t2 + t3*t4 - t3 - t4 + 1, -t2*t3 + t2*t5 + t2 + t3*t4 - t3 - t4*t5 - t4 + t5 + 1]
```
"""
function total_chern_class(
  E::AbstractGKMVectorBundle,
)
  unit = first(gens_cohomRing(baseof(E)))

  return sum(
    chern_classes(E);
    init=zero(unit),
  )
end

function Oscar.chern_class(
  G::AbstractGKMGraph,
  k::Integer,
)
  return chern_class(tangent_bundle(G), k)
end

Oscar.chern_classes(
  G::AbstractGKMGraph,
) = chern_classes(tangent_bundle(G))

total_chern_class(
  G::AbstractGKMGraph,
) = total_chern_class(tangent_bundle(G))

first_chern_class(
  G::AbstractGKMGraph,
) = chern_class(G, 1)
