function _multiplication_matrix_at_q1(G::AbstractGKMGraph, class::GKMClass;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  @req all(e -> chern_number(e, G) > 0, edges(G)) "the graph must be strictly nef"
  class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  # coefficient_ring = parent(first(localize(class).restrictions))
  basis = isnothing(basis) ? [point_class(G, i) for i in 1:num_vertices(G)] : collect(basis)
  n = length(basis)
  n <= num_vertices(G) || throw(DimensionMismatch(
    "the cohomology basis cannot have more than $(num_vertices(G)) elements",
  ))
  all(c -> c isa GKMClass && c.graph === G, basis) || throw(ArgumentError(
    "all basis elements must be GKM classes on the given graph",
  ))

  # If P[j, k] = integral(T_j T_k), first collect
  # S[j, k] = sum_beta GW(class, T_j, T_k). The multiplication matrix is
  # then the solution of P * A = transpose(S), avoiding an explicit inverse.
  # pairing_2 = matrix(coefficient_ring, n, n, [
  #   integrate(basis[j] * basis[k]) for j in 1:n for k in 1:n
  # ])
  # pairing = fast_mode ? _specialize_matrix_at_origin(G, pairing_2) : pairing_2
  pairing = intersection_matrix(basis; equivariant = !fast_mode)
  # invariant_sums = zero_matrix(coefficient_ring, n, n)
  invariant_sums = twisted_intersection_matrix(basis, class; equivariant = !fast_mode)
  symmetric_indices = get_symmetric_indices(n)
  for degree in 0:(2 * valency(G))
    for beta in _effective_classes_with_functional_value(
        GKM_second_homology(G), first_chern_class(G), degree)
      show_progress && println("Computing quantum product in curve class $beta:")
      if iszero(beta)
        # for j in 1:n, k in j:n
        #   invariant = integrate(class * basis[j] * basis[k])
        #   invariant_sums[j, k] += invariant
        #   j == k || (invariant_sums[k, j] += invariant)
        # end
        continue
      end
      # symmetric_indices = [(j, k) for j in 1:n for k in j:n]
      marked_insertions = [
        ev(1, class) * ev(2, basis[j]) * ev(3, basis[k])
        for (j, k) in symmetric_indices
      ]
      invariants = gromov_witten(
        G, beta, 3, marked_insertions; show_bar=show_progress, fast_mode,
      )
      # marked_insertions = [
      #   GKMClass[class, basis[j], basis[k]]
      #   for (j, k) in symmetric_indices
      # ]
      # invariants = gromov_witten_nomarks(
      #   G, beta, marked_insertions; show_bar=show_progress, fast_mode,
      # )
      for ((j, k), invariant) in zip(symmetric_indices, invariants)
        invariant_sums[j, k] += invariant
        j == k || (invariant_sums[k, j] += invariant)
      end
    end
  end
  result = _MP_inv(pairing) * invariant_sums
  # fast_mode ? _specialize_matrix_at_origin(G, result) : result
  return result
end


@doc raw"""
    quantum_product_at_q1(G::AbstractGKMGraph, class;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the matrix of equivariant quantum multiplication by `class` after
setting ``q=1``. The matrix represents

```math
  a \longmapsto (\mathrm{class} \ast a)|_{q=1}
```

in the chosen basis of ``H_T^*(X;\mathbb{Q})``, where `G` is the GKM graph of
``X``. The class can be supplied as a [`GKMClass`](@ref) or as a vector
containing one fixed-point restriction per vertex.

By default, the matrix is expressed in the basis
`[point_class(G, v) for v in 1:num_vertices(G)]`. Pass a collection of
`GKMClass` objects through `basis` to use a different cohomology basis. The
``j``-th column is the product of `class` with the ``j``-th basis element.

The graph must be strictly nef. Otherwise infinitely many effective curve
classes may contribute, so specialization to ``q=1`` need not be defined.
Set `show_progress=true` to display the underlying Gromov--Witten calculations.
Set `fast_mode=true` to use the faster graph-enumeration mode in positive
curve degrees and return the matrix after setting all equivariant parameters
to zero. The chosen basis must remain a basis under this specialization.

# Example

```jldoctest quantum_product_at_q1
julia> P1 = projective_space(GKMGraph, 1);

julia> quantum_product_at_q1(P1, point_class(P1, 1))
[(t1^2 - 2*t1*t2 + t2^2 + 1)//(t1 - t2)    1//(t1 - t2)]
[                         -1//(t1 - t2)   -1//(t1 - t2)]

julia> t1, t2 = gens_coeffRing(P1);

julia> quantum_product_at_q1(P1, [t1, t2])
[(t1^2 - t1*t2 + 1)//(t1 - t2)                    1//(t1 - t2)]
[                -1//(t1 - t2)   (t1*t2 - t2^2 - 1)//(t1 - t2)]
```
"""
function quantum_product_at_q1(G::AbstractGKMGraph, class;
    basis=nothing, show_progress::Bool=true, fast_mode::Bool=true)
  if class isa AbstractVector
    isempty(class) && throw(ArgumentError("class must not be empty"))
    class = localized_class(G, class)
  elseif class isa GKMClass
    class.graph === G || throw(ArgumentError("class belongs to a different GKM graph"))
  else
    throw(ArgumentError("class must be a GKMClass or a vector of fixed-point restrictions"))
  end
  _multiplication_matrix_at_q1(G, class; basis, show_progress, fast_mode)
end


@doc raw"""
    c1_at_q1(G::AbstractGKMGraph;
        basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)

Return the matrix of equivariant quantum multiplication by ``c_1(T_X)`` at
``q=1``. By default, the point classes form the basis; pass a collection of
`GKMClass` objects through `basis` to use another cohomology basis.

The graph must be strictly nef so that only finitely many effective curve
classes contribute. Set `show_progress=true` to show the underlying
Gromov--Witten calculations. Set `fast_mode=true` to use the faster
graph-enumeration mode and return the matrix after setting all equivariant
parameters to zero. The chosen basis must have a well-defined non-equivariant
specialization.

# Example
```jldoctest c1_at_q1
julia> P1 = projective_space(GKMGraph, 1);

julia> c1_at_q1(P1)
[(t1^2 - 2*t1*t2 + t2^2 + 2)//(t1 - t2)                              2//(t1 - t2)]
[                         -2//(t1 - t2)   (-t1^2 + 2*t1*t2 - t2^2 - 2)//(t1 - t2)]
```
"""
function c1_at_q1(G::AbstractGKMGraph;
    basis=nothing, show_progress::Bool=false, fast_mode::Bool=false)
  quantum_product_at_q1(
    G, first_chern_class(G); basis, show_progress, fast_mode,
  )
end
