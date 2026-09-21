"""
    matrix_small_quantum_product(G, class; basis=nothing, beta=nothing, show_bar=true, equivariant=false)
    matrix_small_quantum_product(G, classes::Vector; kwargs...)

Return the matrix of small quantum multiplication by `class`, or a vector of
matrices in the same order as `classes`. With `beta=nothing`, sum the curve
contributions at q=1; otherwise compute only the contribution of `beta`.
The vector overload shares the pairing and batches insertions for each curve.
"""
function matrix_small_quantum_product(G::AbstractGKMGraph, class::GKMClass;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

    return only(matrix_small_quantum_product(G, [class];
    basis = basis, beta = beta, show_bar = show_bar, equivariant = equivariant))
end

function matrix_small_quantum_product(G::AbstractGKMGraph, classes::Vector{GKMClass{GT,R}};
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false) where {GT,R}

  @req all(e -> chern_number(e, G) > 0, edges(G)) "the graph must be strictly nef"
  all(class -> class.graph === G, classes) || throw(ArgumentError("all classes must belong to the given GKM graph"))
  isempty(classes) && return []
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the given graph",
  ))

  inverse_pairing = _MP_inv(intersection_matrix(basis; equivariant = equivariant))

  if !isnothing(beta) && iszero(beta)
    return [inverse_pairing * twisted_intersection_matrix(basis, class; equivariant = equivariant)
            for class in classes]
  end

  symmetric_indices = get_symmetric_indices(n)

  if isnothing(beta)
    invariant_sums = [twisted_intersection_matrix(basis, class; equivariant = equivariant)
                      for class in classes] # degree-zero contributions
    for degree in 0:(2 * valency(G))
      for curve in _effective_classes_with_functional_value(
          GKM_second_homology(G), first_chern_class(G), degree)
        if iszero(curve)
          continue
        end
        show_bar && println("Computing quantum product in curve class $curve:")
        invariants = _invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant)

        for ((i, (j, k)), invariant) in zip(
            ((i, pair) for i in eachindex(classes) for pair in symmetric_indices), invariants)
          invariant_sums[i][j, k] += invariant
          j == k || (invariant_sums[i][k, j] += invariant)
        end
      end
    end
  else
    invariants = _invariants(classes, basis, symmetric_indices, G, beta, show_bar, equivariant)
    Z = zero_matrix(
      equivariant ? parent(first(localize(first(classes)).restrictions)) : QQ, n, n,
    )
    invariant_sums = [deepcopy(Z) for class in classes]

    for ((i, (j, k)), invariant) in zip(
            ((i, pair) for i in eachindex(classes) for pair in symmetric_indices), invariants)
      invariant_sums[i][j, k] += invariant
      j == k || (invariant_sums[i][k, j] += invariant)
    end
  end

  return [inverse_pairing * M for M in invariant_sums]
end

function matrix_small_quantum_product_tangent_class(G::AbstractGKMGraph;
    basis = nothing, beta = nothing, show_bar::Bool = true, equivariant::Bool = false)

  class = first_chern_class(G)
  return matrix_small_quantum_product(G, class;
    basis = basis, beta = beta, show_bar = show_bar, equivariant = equivariant)
end

function _invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant)
  class_products = [
    GKMClass[class, basis[j], basis[k]]
    for class in classes for (j, k) in symmetric_indices
  ]
  return gromov_witten_nomarks(
    G, curve, class_products; show_bar, fast_mode = !equivariant,
  )
end

# function _invariants(classes, basis, symmetric_indices, G, curve, show_bar, equivariant)

#   fast_mode = !equivariant
#   marked_insertions = [
#   ev(1, class) * ev(2, basis[j]) * ev(3, basis[k])
#   for class in classes for (j, k) in symmetric_indices
#   ]
#   invariants = gromov_witten(
#     G, curve, 3, marked_insertions; show_bar = show_bar, fast_mode,
#   )
#   return invariants
# end
