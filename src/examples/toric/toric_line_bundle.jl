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
    copy(labels(G)),
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
  fixed_ray_indices::Vector{Int},
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

"""
    weighted_projective_line_bundle(X::WeightedProjectiveSpace, k::Integer; small_torus=false)

Return the orbifold GKM line bundle `O(k)` on the weighted projective
space `X`.

The default uses the coordinate torus. At the fixed point corresponding to
the `i`-th Cox coordinate, the rational torus fibre weight is
`(k / X.w[i]) e_i`; the local stabilizer acts on the fibre by the character
`k` modulo the local isotropy factors. Rational characters are necessary for
unequal weights and are compatible with the orbifold GKM edge weights.
"""
function weighted_projective_line_bundle(
  X::WeightedProjectiveSpace,
  k::Integer;
  small_torus::Bool = false,
)
  @req !small_torus "O(k) for weighted projective spaces is currently implemented for the coordinate torus; use small_torus=false"

  G = _rationalize_orbifold_graph(
    gkm_graph_of_orbifold_toric(X; small_torus = small_torus),
  )
  M = lattice(G)
  GMtoM = hom(M, M, gens(M))
  fixed_ray_indices = _fixed_ray_indices(X)
  weights = [
    (QQ(k) / X.w[i]) * gens(M)[i]
    for i in fixed_ray_indices
  ]
  fiber_reps = _wps_line_bundle_fiber_representations(G, fixed_ray_indices, k)
  return orbifold_line_bundle(G, M, GMtoM, weights; fiber_reps = fiber_reps)
end

"""
    line_bundle_O(X::WeightedProjectiveSpace, k::Integer; small_torus=false)

Alias for `weighted_projective_line_bundle(X, k; small_torus)`.
"""
function line_bundle_O(
  X::WeightedProjectiveSpace,
  k::Integer;
  small_torus::Bool = false,
)
  return weighted_projective_line_bundle(X, k; small_torus = small_torus)
end
