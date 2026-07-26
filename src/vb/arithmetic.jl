function _weight_element_type(
  ::AbstractGKMVectorBundle{R},
) where {R}
  return AbstractAlgebra.Generic.FreeModuleElem{R}
end

function _check_same_bundle_context(
  E::AbstractGKMVectorBundle,
  F::AbstractGKMVectorBundle,
)
  if !(baseof(E) === baseof(F))
    @warn """
    Vector bundles could be defined on different but isomorphic GKM bases;
    using the base of the first bundle.
    """
  end

  @req E.M == F.M """
  Vector bundles need to have the same character lattice
  """

  @req E.GMtoM == F.GMtoM """
  Vector bundles need to have the same GMtoM map
  """

  return nothing
end

function _check_same_bundle_context(bundles::Tuple)
  @req !isempty(bundles) "Need at least one vector bundle"
  first_bundle = first(bundles)
  for E in Base.tail(bundles)
    _check_same_bundle_context(first_bundle, E)
  end
  return nothing
end

function _mod_character_entry(x, isotropy::Vector{Int}, row::Int)
  return row <= length(isotropy) ? mod(Int(x), isotropy[row]) : Int(x)
end

function _copy_rep_column!(target::ZZMatrix, target_col::Int, source::ZZMatrix, source_col::Int)
  for row in 1:size(target, 1)
    target[row, target_col] = source[row, source_col]
  end
  return nothing
end

function _sum_rep_columns(rep::ZZMatrix, indices, isotropy::Vector{Int})
  values = zeros(Int, size(rep, 1))
  for row in 1:size(rep, 1)
    total = sum(i -> Int(rep[row, i]), indices; init = 0)
    values[row] = _mod_character_entry(total, isotropy, row)
  end
  return values
end

function _set_rep_column!(target::ZZMatrix, target_col::Int, values::Vector{Int})
  for row in 1:size(target, 1)
    target[row, target_col] = values[row]
  end
  return nothing
end

function _zero_line_bundle_like(E::OrbifoldGKMVectorBundle)
  G = baseof(E)
  weights = Matrix{_weight_element_type(E)}(undef, num_vertices(G), 1)
  for v in 1:num_vertices(G)
    weights[v, 1] = zero(E.M)
  end
  fiber_reps = _trivial_fiber_representations(G, 1)
  return orbifold_vector_bundle(G, E.M, E.GMtoM, weights; fiber_reps = fiber_reps)
end

"""
    direct_sum(E::OrbifoldGKMVectorBundle...) -> OrbifoldGKMVectorBundle

Return the direct sum of orbifold GKM vector bundles.
"""
function direct_sum(bundles::OrbifoldGKMVectorBundle...)
  _check_same_bundle_context(bundles)

  first_bundle = first(bundles)
  G = baseof(first_bundle)
  total_rank = sum(rank, bundles; init = 0)
  weights = Matrix{_weight_element_type(first_bundle)}(
    undef,
    num_vertices(G),
    total_rank,
  )
  fiber_reps = [zero_matrix(ZZ, length(G.vertex_isotropy[v].isotropy_group), total_rank) for v in 1:num_vertices(G)]

  offset = 0
  for E in bundles
    for v in 1:num_vertices(G)
      for i in 1:rank(E)
        weights[v, offset + i] = E.weights[v, i]
        _copy_rep_column!(fiber_reps[v], offset + i, E.fiber_reps[v], i)
      end
    end
    offset += rank(E)
  end

  return orbifold_vector_bundle(G, first_bundle.M, first_bundle.GMtoM, weights; fiber_reps = fiber_reps)
end

+(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle) = direct_sum(E, F)

"""
    tensor_product(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle) -> OrbifoldGKMVectorBundle

Return the tensor product of orbifold GKM vector bundles.
"""
function tensor_product(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle)
  _check_same_bundle_context(E, F)

  G = baseof(E)
  tensor_rank = rank(E) * rank(F)
  weights = Matrix{_weight_element_type(E)}(undef, num_vertices(G), tensor_rank)
  fiber_reps = [zero_matrix(ZZ, length(G.vertex_isotropy[v].isotropy_group), tensor_rank) for v in 1:num_vertices(G)]

  for v in 1:num_vertices(G)
    isotropy = G.vertex_isotropy[v].isotropy_group
    col = 1
    for i in 1:rank(E)
      for j in 1:rank(F)
        weights[v, col] = E.weights[v, i] + F.weights[v, j]
        values = zeros(Int, size(fiber_reps[v], 1))
        for row in 1:size(fiber_reps[v], 1)
          total = Int(E.fiber_reps[v][row, i]) + Int(F.fiber_reps[v][row, j])
          values[row] = _mod_character_entry(total, isotropy, row)
        end
        _set_rep_column!(fiber_reps[v], col, values)
        col += 1
      end
    end
  end

  return orbifold_vector_bundle(G, E.M, E.GMtoM, weights; fiber_reps = fiber_reps)
end

*(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle) = tensor_product(E, F)

function _wedge_or_symmetric_product(E::OrbifoldGKMVectorBundle, n::Int, wedged::Bool)
  @req n >= 0 "n must be non-negative"
  @req !wedged || n <= rank(E) "n is greater than the rank"

  n == 0 && return _zero_line_bundle_like(E)
  n == 1 && return E

  index_set = collect(1:rank(E))
  indices = wedged ? collect(powerset(index_set, n, n)) : collect(with_replacement_combinations(index_set, n))
  product_rank = length(indices)
  G = baseof(E)
  weights = Matrix{_weight_element_type(E)}(undef, num_vertices(G), product_rank)
  fiber_reps = [zero_matrix(ZZ, length(G.vertex_isotropy[v].isotropy_group), product_rank) for v in 1:num_vertices(G)]

  for v in 1:num_vertices(G)
    isotropy = G.vertex_isotropy[v].isotropy_group
    for col in 1:product_rank
      weights[v, col] = sum(i -> E.weights[v, i], indices[col]; init = zero(E.M))
      _set_rep_column!(fiber_reps[v], col, _sum_rep_columns(E.fiber_reps[v], indices[col], isotropy))
    end
  end

  return orbifold_vector_bundle(G, E.M, E.GMtoM, weights; fiber_reps = fiber_reps)
end

"""
    wedge_product(E::OrbifoldGKMVectorBundle, n::Integer) -> OrbifoldGKMVectorBundle

Return the exterior power of `E`.
"""
function wedge_product(E::OrbifoldGKMVectorBundle, n::Integer)
  return _wedge_or_symmetric_product(E, Int(n), true)
end

"""
    sym_product(E::OrbifoldGKMVectorBundle, n::Integer) -> OrbifoldGKMVectorBundle

Return the symmetric power `Sym^n E`.
"""
function sym_product(E::OrbifoldGKMVectorBundle, n::Integer)
  return _wedge_or_symmetric_product(E, Int(n), false)
end

function Oscar.det(E::OrbifoldGKMVectorBundle)
  return wedge_product(E, rank(E))
end

function _zero_line_bundle_like(
  E::GKMVectorBundle,
)
  G = baseof(E)

  weights = Matrix{
    _weight_element_type(E)
  }(
    undef,
    num_vertices(G),
    1,
  )

  for v in 1:num_vertices(G)
    weights[v, 1] = zero(E.M)
  end

  return vector_bundle(
    G,
    E.M,
    E.GMtoM,
    weights,
  )
end

"""
    direct_sum(E::GKMVectorBundle...) -> GKMVectorBundle

Return the direct sum of smooth GKM vector bundles.
"""
function direct_sum(
  bundles::GKMVectorBundle...,
)
  _check_same_bundle_context(bundles)

  first_bundle = first(bundles)
  G = baseof(first_bundle)
  total_rank = sum(rank, bundles; init=0)

  weights = Matrix{
    _weight_element_type(first_bundle)
  }(
    undef,
    num_vertices(G),
    total_rank,
  )

  offset = 0

  for E in bundles
    for v in 1:num_vertices(G)
      for i in 1:rank(E)
        weights[v, offset + i] = E.weights[v, i]
      end
    end

    offset += rank(E)
  end

  return vector_bundle(
    G,
    first_bundle.M,
    first_bundle.GMtoM,
    weights,
  )
end

Base.:+(
  E::GKMVectorBundle,
  F::GKMVectorBundle,
) = direct_sum(E, F)

"""
    tensor_product(E::GKMVectorBundle, F::GKMVectorBundle)

Return the tensor product of two smooth GKM vector bundles.
"""
function tensor_product(
  E::GKMVectorBundle,
  F::GKMVectorBundle,
)
  _check_same_bundle_context(E, F)

  G = baseof(E)
  tensor_rank = rank(E) * rank(F)

  weights = Matrix{
    _weight_element_type(E)
  }(
    undef,
    num_vertices(G),
    tensor_rank,
  )

  for v in 1:num_vertices(G)
    column = 1

    for i in 1:rank(E)
      for j in 1:rank(F)
        weights[v, column] =
          E.weights[v, i] + F.weights[v, j]

        column += 1
      end
    end
  end

  return vector_bundle(
    G,
    E.M,
    E.GMtoM,
    weights,
  )
end

Base.:*(
  E::GKMVectorBundle,
  F::GKMVectorBundle,
) = tensor_product(E, F)

function _wedge_or_symmetric_product(
  E::GKMVectorBundle,
  n::Int,
  wedged::Bool,
)
  @req n >= 0 "n must be non-negative"
  @req !wedged || n <= rank(E) """
  n is greater than the rank
  """

  n == 0 && return _zero_line_bundle_like(E)
  n == 1 && return E

  index_set = collect(1:rank(E))

  indices = if wedged
    collect(powerset(index_set, n, n))
  else
    collect(
      with_replacement_combinations(index_set, n),
    )
  end

  G = baseof(E)

  weights = Matrix{
    _weight_element_type(E)
  }(
    undef,
    num_vertices(G),
    length(indices),
  )

  for v in 1:num_vertices(G)
    for column in eachindex(indices)
      weights[v, column] = sum(
        i -> E.weights[v, i],
        indices[column];
        init=zero(E.M),
      )
    end
  end

  return vector_bundle(
    G,
    E.M,
    E.GMtoM,
    weights,
  )
end

function wedge_product(
  E::GKMVectorBundle,
  n::Integer,
)
  return _wedge_or_symmetric_product(
    E,
    Int(n),
    true,
  )
end

function sym_product(
  E::GKMVectorBundle,
  n::Integer,
)
  return _wedge_or_symmetric_product(
    E,
    Int(n),
    false,
  )
end

function Oscar.det(E::GKMVectorBundle)
  return wedge_product(E, rank(E))
end

"""
    total_space(E::GKMVectorBundle) -> GKMGraph

Return the GKM graph of the total space of `E`.

The base flags are transported along the character-lattice map of `E`, and
the fibre weights are added as standalone flags.
"""
function Oscar.total_space(E::GKMVectorBundle{R,V,F}) where {R,V,F}
  G = baseof(E)
  total_flags = Vector{Vector{F}}(undef, num_vertices(G))

  for v in vertices(G)
    total_flags[v] = F[
      F(E.GMtoM(flag.weight))
      for flag in flags(G, v)
    ]
    append!(total_flags[v], (F(E.weights[v, i]) for i in 1:rank(E)))
  end

  total_core = GKMCombinatorialData{R,V,F}(
    deepcopy(graph(G)),
    E.M,
    copy(labels(G)),
    total_flags,
    copy(core(G).edge_flags),
  )
  total_connection = build_gkm_connection(
    total_core;
    connection_type="Total space",
  )
  total_cohomology = create_cohomology(rank(E.M), num_vertices(G))

  return GKMGraph{R,V,F}(
    total_core,
    total_connection,
    total_cohomology,
    nothing,
  )
end

function _orbifold_total_space_connection(E::OrbifoldGKMVectorBundle)
  G = baseof(E)
  base_rank = length(flags(G, first(vertices(G))))
  total_transport = Dict{Edge,Vector{Int}}()
  total_coefficients = Dict{Edge,Vector{QQFieldElem}}()

  for base_edge in edges(G), e in (base_edge, reverse(base_edge))
    total_transport[e] = vcat(
      transport(connection(G))[e],
      base_rank .+ transport(connection(E))[e],
    )
    total_coefficients[e] = vcat(
      QQ.(coefficients(connection(G))[e]),
      QQ.(coefficients(connection(E))[e]),
    )
  end

  return Connection{QQFieldElem}(
    total_transport, total_coefficients, "Orbifold total space",
  )
end

function _orbifold_fiber_weight_order(
  E::OrbifoldGKMVectorBundle,
  i::Int,
)
  return foldl(
    lcm,
    (
      Int(denominator(E.weights[v, i][j]))
      for v in vertices(baseof(E)) for j in 1:rank(E.M)
    );
    init=1,
  )
end

"""
    total_space(E::OrbifoldGKMVectorBundle) -> OrbifoldGKMGraph

Return the orbifold GKM graph of the total space of `E`.  Torus weights and
finite-isotropy representations are both included in the tangent data.
"""
function Oscar.total_space(
  E::OrbifoldGKMVectorBundle{R,V,F},
) where {R,V,F}
  G = baseof(E)
  total_flags = Vector{Vector{F}}(undef, num_vertices(G))
  vertex_isotropy = Vector{OrbifoldVertexIsotropy}(undef, num_vertices(G))

  for v in vertices(G)
    total_flags[v] = F[
      F(E.GMtoM(flag.weight), order_of_generic_stabilizer(flag))
      for flag in flags(G, v)
    ]
    append!(
      total_flags[v],
      (let order = _orbifold_fiber_weight_order(E, i)
        F(order * E.weights[v, i], order)
      end
        for i in 1:rank(E)
      ),
    )

    vertex_isotropy[v] = OrbifoldVertexIsotropy(
      copy(G.vertex_isotropy[v].isotropy_group),
      hcat(G.vertex_isotropy[v].tangent_rep, E.fiber_reps[v]),
    )
  end

  total_core = GKMCombinatorialData{R,V,F}(
    deepcopy(graph(G)), E.M, copy(labels(G)), total_flags,
    copy(core(G).edge_flags),
  )
  total_connection = _orbifold_total_space_connection(E)

  return OrbifoldGKMGraph(
    total_core,
    vertex_isotropy,
    deepcopy(G.flag_isotropy),
    total_connection,
  )
end

"""
    dual(E::GKMVectorBundle) -> GKMVectorBundle

Return the dual of a smooth GKM vector bundle.
"""
function Oscar.dual(E::GKMVectorBundle)
  G = baseof(E)

  weights = Matrix{
    _weight_element_type(E)
  }(
    undef,
    num_vertices(G),
    rank(E),
  )

  for v in 1:num_vertices(G)
    for i in 1:rank(E)
      weights[v, i] = -E.weights[v, i]
    end
  end

  return vector_bundle(
    G,
    E.M,
    E.GMtoM,
    weights,
  )
end

"""
    dual(E::OrbifoldGKMVectorBundle) -> OrbifoldGKMVectorBundle

Return the dual of an orbifold GKM vector bundle. Both torus weights and
isotropy characters are dualized.
"""
function Oscar.dual(E::OrbifoldGKMVectorBundle)
  G = baseof(E)

  weights = Matrix{
    _weight_element_type(E)
  }(
    undef,
    num_vertices(G),
    rank(E),
  )

  fiber_reps = Vector{ZZMatrix}(
    undef,
    num_vertices(G),
  )

  for v in 1:num_vertices(G)
    isotropy = G.vertex_isotropy[v].isotropy_group

    fiber_reps[v] = zero_matrix(
      ZZ,
      size(E.fiber_reps[v], 1),
      rank(E),
    )

    for i in 1:rank(E)
      weights[v, i] = -E.weights[v, i]

      for row in 1:size(E.fiber_reps[v], 1)
        fiber_reps[v][row, i] = mod(
          -Int(E.fiber_reps[v][row, i]),
          isotropy[row],
        )
      end
    end
  end

  return orbifold_vector_bundle(
    G,
    E.M,
    E.GMtoM,
    weights;
    fiber_reps=fiber_reps,
  )
end