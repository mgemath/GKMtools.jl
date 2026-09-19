@doc raw"""
    intersection_matrix(basis::Vector{GKMClass}; equivariant::Bool = false)

Brief description.

# Arguments
- `arguments`: description.

# Returns
- description.

# Examples
```jldoctest
julia> example
```
"""
function Oscar.intersection_matrix(basis::Vector{GKMClass{G,R}}; equivariant::Bool = false)  where {G,R}
    @req !isempty(basis) "The basis is empty"
    optional_class = one(first(basis))
    return _intersection_matrix(basis, optional_class, Val(equivariant))
end

@doc raw"""
    twisted_intersection_matrix(basis::Vector{GKMClass{G,R}}, V::GKMVectorBundle; equivariant::Bool = false) where {G,R}

Brief description.

# Arguments
- `arguments`: description.

# Returns
- description.

# Examples
```jldoctest
julia> example
```
"""
function twisted_intersection_matrix(basis::Vector{GKMClass{G,R}}, V::GKMVectorBundle; equivariant::Bool = false) where {G,R}
    @req !isempty(basis) "The basis is empty"
    optional_class = chern_class(V, rank(V))
    return _intersection_matrix(basis, optional_class, Val(equivariant))
end

function twisted_intersection_matrix(basis, class; equivariant::Bool = false) #where {G,R,GG,RR}
    @req !isempty(basis) "The basis is empty"
    return _intersection_matrix(basis, class, Val(equivariant))
end

function _intersection_matrix(basis, optional_class, ::Val{equivariant}) where {equivariant}
    @req !isempty(basis) "The basis is empty"
    n = length(basis)
    # coefficient_ring = equivariant ? parent(first(basis)).coefficient_ring : parent(first(basis)).localized_coefficient_ring
    coefficient_ring = parent(first(basis)).localized_coefficient_ring
    G = first(basis).graph
    half_pairing = zero_matrix(coefficient_ring, n, n)
    for j in 1:n
        for k in j:n
            half_pairing[j, k] = integrate(optional_class * basis[j] * basis[k])
        end
    end
    
    pairing = equivariant ? half_pairing : GKMtools._specialize_matrix_at_origin(G, half_pairing)

    # symmetrize

    for j in 1:n
        for k in 1:(j-1)
            pairing[j, k] = pairing[k, j]
        end
    end

    return pairing 
end

# function _has_finite_number_of_beta_in_the_product(H2::GKM_H2, class::GKMClass,
#     value::Integer; check::Bool=true)
#   ans = true
#   G = class.graph
#   integers = free_module(ZZ, 1)
#   z = first(gens(integers))
#   for (e, edge_index) in H2.edge_to_gen
#     degree = integrate(class, e)
#     if denominator(degree) == 1 
#       println("the class does not have integral degree on $e")
#       ans = false
#       break
#     end

#     if is_constant(numerator(degree)) 
#       println("the degree on $e is not constant")
#       ans = false
#       break
#     end

#     degree = ZZ(constant_coefficient(numerator(degree)))
#     if check && degree > 0 
#       println("the class must be strictly positive on every edge")
#       ans = false
#       break
#     end
#   end
#   return ans
# end

"""
    _effective_classes_with_functional_value(H2, class, value; check=true)
    _effective_classes_with_functional_value(H2, class; check=true)

For a codimension-one class, enumerate effective `beta` with `class ⋅ beta == value`.
For a homogeneous class of higher codimension, `value` is the remaining
codimension: `codim(class) - c₁(TG) ⋅ beta == value`. Values outside
`0:valency(G)` give no curves. Omitting `value` enumerates the full interval
`0 <= codim(class) - c₁(TG) ⋅ beta <= valency(G)`.

With `check=true`, the enumerating divisor must be strictly positive on every
edge (for higher codimension, this divisor is `c₁(TG)`).
"""
function _effective_classes_with_functional_value(H2::GKM_H2, class::GKMClass;
    check::Bool=true)
  codimension = _homogeneous_degree(class)
  @req !isnothing(codimension) && codimension > 1 "the class must be homogeneous of codimension greater than one when value is omitted"
  return Iterators.flatten(
    _effective_classes_with_functional_value(H2, class, value; check=check)
    for value in 0:valency(class.graph)
  )
end

function _effective_classes_with_functional_value(H2::GKM_H2, class::GKMClass,
    value::Integer; check::Bool=true)
  G = class.graph
  codimension = _homogeneous_degree(class)
  if !isnothing(codimension) && codimension > 1
    0 <= value <= valency(G) || return CurveClass[]
    return _effective_classes_with_functional_value(
      H2, first_chern_class(G), codimension - value; check=check)
  end
  integers = free_module(ZZ, 1)
  z = first(gens(integers))
  edge_values = Vector{typeof(z)}(undef, rank(H2.edge_lattice))
  for (e, edge_index) in H2.edge_to_gen
    degree = integrate(class, e)
    @req denominator(degree) == 1 "the class does not have integral degree on $e"
    @req is_constant(numerator(degree)) "the degree on $e is not constant"
    degree = ZZ(constant_coefficient(numerator(degree)))
    check && @req degree > 0 "the class must be strictly positive on every edge"
    edge_values[edge_index] = degree * z
  end
  edge_to_value = ModuleHomomorphism(H2.edge_lattice, integers, edge_values)
  H2_to_value = ModuleHomomorphism(H2.H2, integers,
    [edge_to_value(preimage(H2.quotient, g)) for g in gens(H2.H2)])
  success, beta0 = has_preimage_with_preimage(H2_to_value, value * z)
  success || return CurveClass[]
  K, inclusion = kernel(H2_to_value)
  inclusion_matrix = transpose(matrix(inclusion))
  dual_rays = rays(H2.dual_cone)
  inequalities = matrix(QQ, length(dual_rays), rank(H2.H2),
    [dual_rays[i][j] for i in eachindex(dual_rays) for j in 1:rank(H2.H2)])
  rhs = inequalities * [beta0[i] for i in 1:rank(H2.H2)]
  polytope = polyhedron(-inequalities * inclusion_matrix, rhs)
  (beta0 + inclusion(K([p[i] for i in 1:rank(K)])) for p in lattice_points(polytope))
end


function _specialize_matrix_at_origin(G::AbstractGKMGraph, M)
  origin = zeros(Int, rank_torus(G))
  matrix(QQ, nrows(M), ncols(M), [begin
    numerator_value = evaluate(numerator(M[i, j]), origin)
    denominator_value = evaluate(denominator(M[i, j]), origin)
    iszero(denominator_value) && throw(ArgumentError(
      "the multiplication matrix has no entrywise non-equivariant limit " *
      "in the chosen basis; use a basis that remains a basis at t=0",
    ))
    QQ(numerator_value) / QQ(denominator_value)
  end for i in 1:nrows(M) for j in 1:ncols(M)])
end

function _MP_inv(A)
  m, n = nrows(A), ncols(A)
  coefficient_ring = base_ring(A)
  r, reduced = rref(A)
  iszero(r) && return zero_matrix(coefficient_ring, n, m)

  # Choose the pivot columns to obtain a rank factorization A = C * F,
  # with C of full column rank and F of full row rank.
  pivots = Vector{Int}(undef, r)
  for i in 1:r
    pivot = findfirst(j -> !iszero(reduced[i, j]), 1:n)
    isnothing(pivot) && error("failed to find a pivot in a nonzero RREF row")
    pivots[i] = pivot
  end
  C = A[:, pivots]
  F = transpose(solve(transpose(C), transpose(A)))

  # For A = C*F, the Moore--Penrose inverse is F^+*C^+.
  Ft = transpose(F)
  Ct = transpose(C)
  F_plus = solve(F * Ft, Ft)
  C_plus = transpose(solve(transpose(Ct * C), C))
  F_plus * C_plus
end

function get_symmetric_indices(n)
  return [(j, k) for j in 1:n for k in j:n]
end