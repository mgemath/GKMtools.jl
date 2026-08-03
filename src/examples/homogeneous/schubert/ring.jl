@doc raw"""
    schubert_cohomology_ring(G; degree_convention=:cohomological)

Construct the ordinary rational cohomology ring of the homogeneous GKM variety
`G = G/P` as an Oscar graded quotient ring, using the polynomial Schubert basis
and Billey's restrictions. 
See [Gradings](https://docs.oscar-system.org/stable/CommutativeAlgebra/rings/#Gradings).

The result is `(H, sigma)`. Here `H` is the graded quotient and
`sigma[d + 1]` is the vector of quotient elements corresponding to the
codimension-`d` Schubert basis. Thus `sigma[1] == [one(H)]`.

By default a codimension-`d` class has degree `2d`. Set
`degree_convention=:codimension` to give it degree `d` instead.

The presentation uses one variable for every positive-codimension Schubert
class and relations for their complete multiplication table. It is therefore
canonical relative to the fixed-point ordering, though not generator-minimal.

The Schubert class corresponding to `g` is denoted by `σ_g`.

# Examples
Let us start with the projective space $\mathbb{P}^n$.
```jldoctest
julia> n = 4;

julia> R = root_system(:A, n);

julia> P4 = generalized_gkm_flag(R, 2:n) # P^4
GKM graph with 5 nodes, valency 4 and axial function:
s1 -> id => (-1, 1, 0, 0, 0)
s2*s1 -> id => (-1, 0, 1, 0, 0)
s2*s1 -> s1 => (0, -1, 1, 0, 0)
s3*s2*s1 -> id => (-1, 0, 0, 1, 0)
s3*s2*s1 -> s1 => (0, -1, 0, 1, 0)
s3*s2*s1 -> s2*s1 => (0, 0, -1, 1, 0)
s4*s3*s2*s1 -> id => (-1, 0, 0, 0, 1)
s4*s3*s2*s1 -> s1 => (0, -1, 0, 0, 1)
s4*s3*s2*s1 -> s2*s1 => (0, 0, -1, 0, 1)
s4*s3*s2*s1 -> s3*s2*s1 => (0, 0, 0, -1, 1)
Birkhoff-Grothendieck connection for GKM graph with 5 nodes and valency 4

julia> (H, sigma) = schubert_cohomology_ring(P4);

julia> H
Quotient
  of multivariate polynomial ring in 4 variables over QQ graded by
    σ_s1 -> [2]
    σ_s2*s1 -> [4]
    σ_s3*s2*s1 -> [6]
    σ_s4*s3*s2*s1 -> [8]
  by ideal with 10 generators

julia> sigma
5-element Vector{Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [1]
 [σ_s1]
 [σ_s2*s1]
 [σ_s3*s2*s1]
 [σ_s4*s3*s2*s1]

julia> σ_s1 = sigma[2][1]; # hyperplane class

julia> foreach(i -> println("(σ_s1)^$i = ", σ_s1^i), 0:5)
(σ_s1)^0 = 1
(σ_s1)^1 = σ_s1
(σ_s1)^2 = σ_s2*s1
(σ_s1)^3 = σ_s3*s2*s1
(σ_s1)^4 = σ_s4*s3*s2*s1
(σ_s1)^5 = 0
```
Let us consider the Grassmannian $G(2, 4)$.
```jldoctest
julia> G24 = generalized_gkm_flag(root_system(:A, 3), [1, 3])
GKM graph with 6 nodes, valency 4 and axial function:
s2 -> id => (0, -1, 1, 0)
s1*s2 -> id => (-1, 0, 1, 0)
s1*s2 -> s2 => (-1, 1, 0, 0)
s3*s2 -> id => (0, -1, 0, 1)
s3*s2 -> s2 => (0, 0, -1, 1)
s1*s3*s2 -> id => (-1, 0, 0, 1)
s1*s3*s2 -> s1*s2 => (0, 0, -1, 1)
s1*s3*s2 -> s3*s2 => (-1, 1, 0, 0)
s2*s1*s3*s2 -> s2 => (-1, 0, 0, 1)
s2*s1*s3*s2 -> s1*s2 => (0, -1, 0, 1)
s2*s1*s3*s2 -> s3*s2 => (-1, 0, 1, 0)
s2*s1*s3*s2 -> s1*s3*s2 => (0, -1, 1, 0)
Birkhoff-Grothendieck connection for GKM graph with 6 nodes and valency 4

julia> (H, sigma) = schubert_cohomology_ring(G24);

julia> H
Quotient
  of multivariate polynomial ring in 5 variables over QQ graded by
    σ_s2 -> [2]
    σ_s1*s2 -> [4]
    σ_s3*s2 -> [4]
    σ_s1*s3*s2 -> [6]
    σ_s2*s1*s3*s2 -> [8]
  by ideal with 15 generators
```
Those are all Schubert classes of the cohomology ring of $G(2, 4)$, where `σ_s2*s1*s3*s2` is the point class and `σ_s2` is the Plücker class.
"""
function schubert_cohomology_ring(
  G::AbstractGKMGraph{R,V,F};
  degree_convention::Symbol=:cohomological,
) where {R,V<:GeneralizedFlagVertex,F}
  degree_factor = if degree_convention == :cohomological
    2
  elseif degree_convention == :codimension
    1
  else
    throw(ArgumentError(
      "degree_convention must be :cohomological or :codimension",
    ))
  end

  vertex_order = sort(
    collect(1:num_vertices(G));
    by=v -> (length(flag(vertices_structure(G)[v])), v),
  )
  codimensions = [length(flag(vertices_structure(G)[v])) for v in vertex_order]
  polynomial_basis = [
    schubert_basis(G, v; representation=:polynomial)
    for v in vertex_order
  ]

  # The identity is the unique codimension-zero Schubert class. Every other
  # basis element is used as a presentation variable.
  codimensions[1] == 0 || error("the Schubert ordering has no identity class")
  count(==(0), codimensions) == 1 || error("the Schubert basis has multiple degree-zero classes")
  positive_positions = collect(2:length(polynomial_basis))
  variable_names = [
    _schubert_variable_name(label(G, vertex_order[k]))
    for k in positive_positions
  ]
  variable_degrees = degree_factor .* codimensions[positive_positions]
  uses_dummy_variable = isempty(positive_positions)
  if uses_dummy_variable
    variable_names = ["_zero"]
    variable_degrees = [degree_factor]
  end

  raw_ring, _ = polynomial_ring(QQ, variable_names)
  graded_ring, variables = grade(raw_ring, variable_degrees)
  multiplication = _ordinary_schubert_multiplication(
    G,
    polynomial_basis,
    vertex_order,
    codimensions,
  )

  relations = uses_dummy_variable ? [variables[1]] : elem_type(graded_ring)[]
  for i in positive_positions
    xi = variables[i - 1]
    for j in i:length(polynomial_basis)
      xj = variables[j - 1]
      rhs = zero(graded_ring)
      for k in eachindex(polynomial_basis)
        coefficient = multiplication[i, j][k]
        iszero(coefficient) && continue
        rhs += coefficient * (k == 1 ? one(graded_ring) : variables[k - 1])
      end
      push!(relations, xi * xj - rhs)
    end
  end

  quotient, quotient_map = quo(graded_ring, ideal(graded_ring, relations))
  flat_basis = vcat(
    [one(quotient)],
    [quotient_map(variables[k - 1]) for k in positive_positions],
  )
  grouped_basis = [elem_type(quotient)[] for _ in 0:maximum(codimensions)]
  for k in eachindex(flat_basis)
    push!(grouped_basis[codimensions[k] + 1], flat_basis[k])
  end

  return quotient, grouped_basis
end

function schubert_cohomology_ring(::AbstractGKMGraph; degree_convention::Symbol=:cohomological)
  degree_convention in (:cohomological, :codimension) || throw(ArgumentError(
    "degree_convention must be :cohomological or :codimension",
  ))
  throw(ArgumentError(
    "a Schubert cohomology presentation requires Weyl-backed generalized flag vertices",
  ))
end

function _schubert_variable_name(vertex_label::String)
  # cleaned = replace(vertex_label, r"[^A-Za-z0-9]+" => "_")
  # return "σ_" * strip(cleaned, '_')
  return "σ_" * vertex_label
end

# Compute all ordinary structure constants. Products are first expanded over
# S = H_T^*(pt) using triangular Billey restrictions, then specialized at t=0.
function _ordinary_schubert_multiplication(
  G::AbstractGKMGraph,
  basis::Vector{<:GKMClass},
  vertex_order::Vector{Int},
  codimensions::Vector{Int},
)
  n = length(basis)
  table = Matrix{Vector{QQFieldElem}}(undef, n, n)
  zeros_at_origin = fill(QQ(0), ngens(equivariant_coefficient_ring(G)))

  for i in 1:n
    for j in i:n
      equivariant_coefficients = _expand_in_schubert_basis(
        basis[i] * basis[j],
        basis,
        vertex_order,
      )
      ordinary_coefficients = QQFieldElem[
        QQ(evaluate(c, zeros_at_origin))
        for c in equivariant_coefficients
      ]

      expected_codimension = codimensions[i] + codimensions[j]
      for k in 1:n
        if codimensions[k] != expected_codimension && !iszero(ordinary_coefficients[k])
          error("ordinary Schubert multiplication did not preserve degree")
        end
      end

      table[i, j] = ordinary_coefficients
      table[j, i] = ordinary_coefficients
    end
  end
  return table
end

function _expand_in_schubert_basis(
  product::GKMClass,
  basis::Vector{<:GKMClass},
  vertex_order::Vector{Int},
)
  S = coefficient_ring(parent(product))
  coefficients = fill(zero(S), length(basis))

  for k in eachindex(basis)
    vertex = vertex_order[k]
    residual = product[vertex]
    for j in 1:k-1
      iszero(coefficients[j]) && continue
      residual -= coefficients[j] * basis[j][vertex]
    end

    diagonal = basis[k][vertex]
    success, quotient = divides(residual, diagonal)
    success || error(
      "Schubert product has a non-polynomial coefficient at $(label(graph(parent(product)), vertex))",
    )
    coefficients[k] = quotient
  end

  # This also detects a failure of the triangular solve away from the diagonal.
  reconstructed = zero(parent(product))
  for k in eachindex(basis)
    iszero(coefficients[k]) || (reconstructed += coefficients[k] * basis[k])
  end
  reconstructed == product || error("Schubert restrictions did not span the product")
  return coefficients
end
