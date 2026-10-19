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
    """ _file=""
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

# Examples
```jldoctest
julia> G = projective_space(GKMGraph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = hom(lattice(G), M, [g[1], g[2], g[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]]);

julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = direct_sum(V1, V2)
GKM vector bundle of rank 2 over 3 vertices with weights:
1: (0, 0, 0, 1), (1, 0, 0, 0)
2: (0, 0, 0, 1), (0, 1, 0, 0)
3: (0, 0, 0, 1), (0, 0, 1, 0)
```
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

# """
#     tensor_product(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle) -> OrbifoldGKMVectorBundle

# Return the tensor product of orbifold GKM vector bundles.
# """
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

# Examples
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKMGraph, 2, 5);

julia> wedge_product(Q, 2)
GKM vector bundle of rank 3 over 10 vertices with weights:
12: (0, 0, -1, -1, 0), (0, 0, -1, 0, -1), (0, 0, 0, -1, -1)
13: (0, -1, 0, -1, 0), (0, -1, 0, 0, -1), (0, 0, 0, -1, -1)
14: (0, -1, -1, 0, 0), (0, -1, 0, 0, -1), (0, 0, -1, 0, -1)
15: (0, -1, -1, 0, 0), (0, -1, 0, -1, 0), (0, 0, -1, -1, 0)
23: (-1, 0, 0, -1, 0), (-1, 0, 0, 0, -1), (0, 0, 0, -1, -1)
24: (-1, 0, -1, 0, 0), (-1, 0, 0, 0, -1), (0, 0, -1, 0, -1)
25: (-1, 0, -1, 0, 0), (-1, 0, 0, -1, 0), (0, 0, -1, -1, 0)
34: (-1, -1, 0, 0, 0), (-1, 0, 0, 0, -1), (0, -1, 0, 0, -1)
35: (-1, -1, 0, 0, 0), (-1, 0, 0, -1, 0), (0, -1, 0, -1, 0)
45: (-1, -1, 0, 0, 0), (-1, 0, -1, 0, 0), (0, -1, -1, 0, 0)
```
"""
function wedge_product(E::OrbifoldGKMVectorBundle, n::Integer)
  return _wedge_or_symmetric_product(E, Int(n), true)
end

"""
    sym_product(E::OrbifoldGKMVectorBundle, n::Integer) -> OrbifoldGKMVectorBundle

Return the symmetric power `Sym^n E`.

# Examples
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKMGraph, 3, 5);

julia> sym_product(Q, 2)
GKM vector bundle of rank 3 over 10 vertices with weights:
123: (0, 0, 0, -2, 0), (0, 0, 0, -1, -1), (0, 0, 0, 0, -2)
124: (0, 0, -2, 0, 0), (0, 0, -1, 0, -1), (0, 0, 0, 0, -2)
125: (0, 0, -2, 0, 0), (0, 0, -1, -1, 0), (0, 0, 0, -2, 0)
134: (0, -2, 0, 0, 0), (0, -1, 0, 0, -1), (0, 0, 0, 0, -2)
135: (0, -2, 0, 0, 0), (0, -1, 0, -1, 0), (0, 0, 0, -2, 0)
145: (0, -2, 0, 0, 0), (0, -1, -1, 0, 0), (0, 0, -2, 0, 0)
234: (-2, 0, 0, 0, 0), (-1, 0, 0, 0, -1), (0, 0, 0, 0, -2)
235: (-2, 0, 0, 0, 0), (-1, 0, 0, -1, 0), (0, 0, 0, -2, 0)
245: (-2, 0, 0, 0, 0), (-1, 0, -1, 0, 0), (0, 0, -2, 0, 0)
345: (-2, 0, 0, 0, 0), (-1, -1, 0, 0, 0), (0, -2, 0, 0, 0)
```
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

# Examples
```jldoctest
julia> G = projective_space(GKMGraph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = hom(lattice(G), M, [g[1], g[2], g[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]]);

julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = direct_sum(V1, V2)
GKM vector bundle of rank 2 over 3 vertices with weights:
1: (0, 0, 0, 1), (1, 0, 0, 0)
2: (0, 0, 0, 1), (0, 1, 0, 0)
3: (0, 0, 0, 1), (0, 0, 1, 0)
```
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
```jldoctest
julia> E = line_bundle_O(3, 1)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (1, 0, 0, 0)
2: (0, 1, 0, 0)
3: (0, 0, 1, 0)
4: (0, 0, 0, 1)

julia> F = line_bundle_O(3, 3)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (1, 0, 0, 0)
2: (-2, 3, 0, 0)
3: (-2, 0, 3, 0)
4: (-2, 0, 0, 3)

julia> tensor_product(E, F)
┌ Warning: Vector bundles could be defined on different but isomorphic GKM bases;
│ using the base of the first bundle.
└ @ GKMtools :12
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (2, 0, 0, 0)
2: (-2, 4, 0, 0)
3: (-2, 0, 4, 0)
4: (-2, 0, 0, 4)

julia> E*F
┌ Warning: Vector bundles could be defined on different but isomorphic GKM bases;
│ using the base of the first bundle.
└ @ GKMtools :12
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (2, 0, 0, 0)
2: (-2, 4, 0, 0)
3: (-2, 0, 4, 0)
4: (-2, 0, 0, 4)
```
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

# Examples
```jldoctest
julia> G = projective_space(GKMGraph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Birkhoff-Grothendieck connection for GKM graph with 3 nodes and valency 2

julia> M = free_module(ZZ, 4);

julia> GMtoM = hom(lattice(G), M, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]])
GKM vector bundle of rank 1 over 3 vertices with weights:
1: (0, 0, 0, 1)
2: (0, 0, 0, 1)
3: (0, 0, 0, 1)

julia> T = total_space(V)
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
Standalone flags:
Total space connection for GKM graph with 3 nodes and valency 3

julia> valency(T)
3

julia> is_compact(T)
false
```

```jldoctest
julia> F = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)
Algorithmic connection for GKM graph with 6 nodes and valency 3

julia> X = total_space(dual(tangent_bundle(F)));

julia> print_curve_classes(X)
2 -> 1: (0, 1), Chern number: 0
3 -> 2: (-1, 1), Chern number: 0
4 -> 1: (1, 0), Chern number: 0
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 0
5 -> 4: (-1, 1), Chern number: 0
6 -> 1: (1, 1), Chern number: 0
6 -> 3: (1, 0), Chern number: 0
6 -> 5: (0, 1), Chern number: 0
```
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
    copy(vertices_structure(G)),
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

# Examples
```jldoctest
julia> G = projective_space(GKMGraph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)
Birkhoff-Grothendieck connection for GKM graph with 3 nodes and valency 2

julia> M = free_module(ZZ, 4);

julia> GMtoM = hom(lattice(G), M, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]])
GKM vector bundle of rank 1 over 3 vertices with weights:
1: (0, 0, 0, 1)
2: (0, 0, 0, 1)
3: (0, 0, 0, 1)

julia> T = total_space(V)
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
Standalone flags:
Total space connection for GKM graph with 3 nodes and valency 3

julia> valency(T)
3

julia> is_compact(T)
false
```

```jldoctest
julia> F = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)
Algorithmic connection for GKM graph with 6 nodes and valency 3

julia> X = total_space(dual(tangent_bundle(F)));

julia> print_curve_classes(X)
2 -> 1: (0, 1), Chern number: 0
3 -> 2: (-1, 1), Chern number: 0
4 -> 1: (1, 0), Chern number: 0
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 0
5 -> 4: (-1, 1), Chern number: 0
6 -> 1: (1, 1), Chern number: 0
6 -> 3: (1, 0), Chern number: 0
6 -> 5: (0, 1), Chern number: 0
```
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
    deepcopy(graph(G)), E.M, copy(vertices_structure(G)), total_flags,
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

# Examples
```jldoctest dual_vector_bundles
julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5));

julia> M = free_module(ZZ, 5);

julia> g = gens(M);

julia> GMtoM = hom(lattice(G), M, [g[1], g[2], g[3], g[4]]);

julia> L = line_bundle(G, M, GMtoM, [g[5], g[5], g[5], g[5]])
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (0, 0, 0, 0, 1)
2: (0, 0, 0, 0, 1)
3: (0, 0, 0, 0, 1)
4: (0, 0, 0, 0, 1)

julia> dual(L)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (0, 0, 0, 0, -1)
2: (0, 0, 0, 0, -1)
3: (0, 0, 0, 0, -1)
4: (0, 0, 0, 0, -1)
```
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

# Examples
```jldoctest dual_vector_bundles
julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5));

julia> M = free_module(ZZ, 5);

julia> g = gens(M);

julia> GMtoM = hom(lattice(G), M, [g[1], g[2], g[3], g[4]]);

julia> L = line_bundle(G, M, GMtoM, [g[5], g[5], g[5], g[5]])
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (0, 0, 0, 0, 1)
2: (0, 0, 0, 0, 1)
3: (0, 0, 0, 0, 1)
4: (0, 0, 0, 0, 1)

julia> dual(L)
GKM vector bundle of rank 1 over 4 vertices with weights:
1: (0, 0, 0, 0, -1)
2: (0, 0, 0, 0, -1)
3: (0, 0, 0, 0, -1)
4: (0, 0, 0, 0, -1)
```
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

function inv(V::AbstractGKMVectorBundle)
  return dual(V)
end

function ^(V::AbstractGKMVectorBundle, n::Integer)
  return sym_product(V, n)
end

@doc raw"""
    projectivization(E::GKMVectorBundle) -> GKMGraph

Return the GKM graph of the projectivization of `E`. Its vertices are the
one-dimensional weight spaces in the fibres of `E`; the vertex corresponding
to the `i`-th weight over a base vertex `v` is labelled `[v]_i`.

# Example
```jldoctest
julia> E = tangent_bundle(projective_space(GKMGraph, 1));

julia> P = projectivization(E);

julia> (num_vertices(P), valency(P))
(2, 1)

julia> [label(P, v) for v in vertices(P)]
2-element Vector{String}:
 "[1]_1"
 "[2]_1"
```
"""
function Oscar.projectivization(E::GKMVectorBundle)
  G = baseof(E)
  total = total_space(E)
  base_flags = [collect(1:length(flags(G, v))) for v in vertices(G)]
  zero_section = _subgraph_from_selected_flags(total, collect(vertices(total)), base_flags)
  P = subgraph(blowup(zero_section))
  projective_labels = vertices_structure(P)
  k = 0
  for v in vertices(G), i in 1:rank(E)
    k += 1
    old = projective_labels[k]
    projective_labels[k] = BlowupVertex("[$(label(G, v))]_$i", old.old_vertex)
  end
  return P
end


"""
    *(pair::Tuple{GKMGraph{R}, GKMVectorBundle},
      pairs::Tuple{GKMGraph{R}, GKMVectorBundle}...) where {R}

Return `(G1 * G2 * ..., (v1, v2, ...))` for one or more graph–bundle
pairs, where `vi` is the pullback of the i-th bundle to the product graph.
Equivalently, each `vi` is the external tensor product with the trivial
line bundles on all other factors. The pullbacks share the direct sum of
all bundle character lattices. For one pair, return `(G1, (V1,))`.
"""
function Base.:*(
  pair::Tuple{GKMGraph{R}, GKMVectorBundle},
  pairs::Tuple{GKMGraph{R}, GKMVectorBundle}...,
) where {R}
  factors = (pair, pairs...)
  for (G, V) in factors
    @req baseof(V) === G "Each vector bundle must be defined on its paired graph"
  end
  isempty(pairs) && return (pair[1], (pair[2],))

  G = foldl(*, map(first, factors))
  bundles = map(last, factors)
  M = free_module(base_ring(first(bundles).M), sum(V -> rank(V.M), bundles))
  lattice_offset = 0
  embeddings = map(bundles) do V
    f = hom(V.M, M, gens(M)[lattice_offset + 1:lattice_offset + rank(V.M)])
    lattice_offset += rank(V.M)
    return f
  end
  images = AbstractAlgebra.Generic.FreeModuleElem{R}[]
  for ((Gi, Vi), f) in zip(factors, embeddings)
    append!(images, (f(Vi.GMtoM(g)) for g in gens(lattice(Gi))))
  end
  GMtoM = hom(lattice(G), M, images)

  # In the left-associated product, the first factor's vertex varies fastest.
  stride = 1
  pullbacks = map(factors, embeddings) do (Gi, Vi), f
    weights = Matrix{_weight_element_type(Vi)}(undef, num_vertices(G), rank(Vi))
    for v in vertices(G)
      vi = mod(div(v - 1, stride), num_vertices(Gi)) + 1
      for i in 1:rank(Vi)
        weights[v, i] = f(Vi.weights[vi, i])
      end
    end
    stride *= num_vertices(Gi)
    return vector_bundle(G, M, GMtoM, weights)
  end
  return G, pullbacks
end
