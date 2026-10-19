function _maximal_cone_rows(X::WeightedProjectiveSpace)
  d = dim(X)
  return findall(i -> sum(X.incidence_matrix[i, :]) == d, 1:size(X.incidence_matrix, 1))
end

function _fixed_ray_indices(X::WeightedProjectiveSpace)
  return [
    only(findall(j -> !X.incidence_matrix[i, j], 1:n_rays(X)))
    for i in _maximal_cone_rows(X)
  ]
end

function _rationalize_orbifold_graph(
  G::OrbifoldGKMGraph{R,V,OrbifoldToricFlagWeight{R}},
) where {R,V}
  R === QQFieldElem && return G

  M = free_module(QQ, rank(lattice(G)))
  basis = gens(M)
  lift_weight(w) = sum(
    i -> QQ(w[i]) * basis[i],
    eachindex(basis);
    init = zero(M),
  )
  rational_flags = [
    [
      OrbifoldToricFlagWeight{QQFieldElem}(
        lift_weight(flag.weight),
        order_of_generic_stabilizer(flag),
      )
      for flag in flags(G, v)
    ]
    for v in vertices(G)
  ]
  data = GKMCombinatorialData{
    QQFieldElem,
    V,
    OrbifoldToricFlagWeight{QQFieldElem},
  }(
    graph(G),
    M,
    copy(vertices_structure(G)),
    rational_flags,
    copy(core(G).edge_flags),
  )
  con = build_gkm_connection(data; connection_type = connection_type(G))
  return OrbifoldGKMGraph(
    data,
    copy(G.vertex_isotropy),
    copy(G.flag_isotropy),
    con,
  )
end

function _wps_line_bundle_fiber_representations(
  G::AbstractOrbifoldGKMGraph,
  k::Integer,
)
  reps = Vector{ZZMatrix}(undef, num_vertices(G))
  for v in 1:num_vertices(G)
    isotropy = G.vertex_isotropy[v].isotropy_group
    reps[v] = zero_matrix(ZZ, length(isotropy), 1)
    for a in 1:length(isotropy)
      reps[v][a, 1] = mod(k, isotropy[a])
    end
  end
  return reps
end

# @doc raw"""
#     weighted_projective_line_bundle(X::WeightedProjectiveSpace, k::Integer; big_torus=false)

# Return the orbifold GKM line bundle $\mathcal{O}(k)$ on the weighted projective
# space `X`.

# The default uses the coordinate torus. At the fixed point corresponding to
# the `i`-th Cox coordinate, the rational torus fibre weight is
# `-(k / X.w[i]) e_i`; the local stabilizer acts on the fibre by the character
# `k` modulo the local isotropy factors. Rational characters are necessary for
# unequal weights and are compatible with the orbifold GKM edge weights.
# """
function weighted_projective_line_bundle(
  X::WeightedProjectiveSpace,
  k::Integer;
  big_torus::Bool = false,
)
  @req !big_torus "O(k) for weighted projective spaces is currently implemented for the coordinate torus; use big_torus=false"

  G = _rationalize_orbifold_graph(
    gkm_graph_of_orbifold_toric(X; big_torus = big_torus),
  )
  M = lattice(G)
  GMtoM = hom(M, M, gens(M))
  fixed_ray_indices = _fixed_ray_indices(X)
  weights = [
    (QQ(-k) / X.w[i]) * gens(M)[i]
    for i in fixed_ray_indices
  ]
  fiber_reps = _wps_line_bundle_fiber_representations(G, k)
  return orbifold_line_bundle(G, M, GMtoM, weights; fiber_reps = fiber_reps)
end

function _line_bundle_O(
  G::GKMGraph,
  k::Integer,
)
  M = lattice(core(G))
  GMtoM = hom(M, M, gens(M))
  weights = Vector{eltype(gens(M))}(undef, num_vertices(G))
  weights[1] = gens(M)[1]

  for v in 2:num_vertices(G)
    weights[v] = weights[1] - k * weight(G, Edge(1, v))
  end

  return line_bundle(G, M, GMtoM, weights)
end

@doc raw"""
    line_bundle_O(n::Integer, k::Integer; big_torus=false)

Return the smooth GKM line bundle $\mathcal{O}(k)$ on the projective `n`-space.
The parameter `big_torus=false`, the default option, means that the line bundle is constructed with the same torus action on the projective space. 
If `big_torus=false`, the default option, the torus acting is the sum of the torus on the projective space and a torus of dimension equal to the rank, acting on fibers.
```jldoctest
julia> E = line_bundle_O(3, 2, big_torus = true)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (0, 0, 0, 0, 1)
2: (-2, 2, 0, 0, 1)
3: (-2, 0, 2, 0, 1)
4: (-2, 0, 0, 2, 1)

julia> E = line_bundle_O(3, 2)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (1, 0, 0, 0)
2: (-1, 2, 0, 0)
3: (-1, 0, 2, 0)
4: (-1, 0, 0, 2)
```
"""
function line_bundle_O(
  n::Integer,
  k::Integer;
  big_torus::Bool=false,
)
  @req n >= 1 "The dimension must be positive"

  if big_torus
    return vector_bundle_O(n, [k]; big_torus = true)
  end

  return _line_bundle_O(
    projective_space(GKMGraph, Int(n)),
    k,
  )
end

@doc raw"""
    line_bundle_O(n::Integer; big_torus=false)

Return the smooth GKM line bundle $\mathcal{O}(1)$ on the projective `n`-space.
"""
function line_bundle_O(
  n::Integer,
  big_torus::Bool=false,
)
  @req n >= 1 "The dimension must be positive"
  @req !big_torus "line_bundle_O(n, k) currently uses the coordinate torus; use big_torus=false"

  return _line_bundle_O(
    projective_space(GKMGraph, Int(n)),
    1,
  )
end

@doc raw"""
    line_bundle_O(X::WeightedProjectiveSpace, k::Integer; big_torus=false)

Alias for `weighted_projective_line_bundle(X, k; big_torus)`.
"""
function line_bundle_O(
  X::WeightedProjectiveSpace,
  k::Integer;
  big_torus::Bool = false,
)
  return weighted_projective_line_bundle(X, k; big_torus = big_torus)
end

@doc raw"""
    line_bundle_O(X::WeightedProjectiveSpace, divisor::AbstractVector{<:Integer}; big_torus=false)

Return the line bundle associated to a torus-invariant Cox divisor
`sum(divisor[i] * D_i)`.

For a weighted projective space `P(w_1,...,w_n)`, the divisor `D_i`
corresponds to `O(w_i)`, so this is `O(sum(divisor .* X.w))`.
"""
function line_bundle_O(
  X::WeightedProjectiveSpace,
  divisor::AbstractVector{<:Integer};
  big_torus::Bool = false,
)
  @req length(divisor) == length(X.w) """
  Divisor vector must have one coefficient for each Cox ray
  """

  k = sum(i -> Int(divisor[i]) * X.w[i], eachindex(divisor); init = 0)
  return weighted_projective_line_bundle(X, k; big_torus = big_torus)
end
