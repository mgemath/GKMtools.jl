"""
    matrix_small_quantum_product(G, class; basis=nothing, beta=nothing, show_bar=true, equivariant=false)
    matrix_small_quantum_product(G, classes::Vector; kwargs...)

Return the matrix of small quantum multiplication by `class`, or a vector of
matrices in the same order as `classes`. With `beta=nothing`, sum the curve
contributions at q=1; otherwise compute only the contribution of `beta`.
Use `basis` to specify a cohomology basis; otherwise the point classes are used.
Use `show_bar=false` to suppress the progress bar. 
If `equivariant=true`, compute the equivariant quantum product. 
If `equivariant=false`, compute the non-equivariant quantum product, in that case a basis should be provided.

```jldoctest
julia> L = line_bundle_O(1, 1)
GKM vector bundle of rank 1 over 2 vertices with weights:
1: (1, 0)
2: (0, 1)

julia> G = baseof(L)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)
Birkhoff-Grothendieck connection for GKM graph with 2 nodes and valency 1

julia> matrix_small_quantum_product(G, first_chern_class(G); equivariant = true, show_bar = false)
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]

julia> basis = [chern_class(L, i-1) for i in 1:num_vertices(G)];

julia> matrix_small_quantum_product(G, first_chern_class(G); basis = basis, equivariant = false, show_bar = false)
[0   2]
[2   0]

julia> matrix_small_quantum_product(G, first_chern_class(G); basis = basis, equivariant = true, show_bar = false)
[-t1 - t2   -2*t1*t2 + 2]
[       2        t1 + t2]
```
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
