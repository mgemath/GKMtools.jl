@doc raw"""
    small_quantum_cohomology_ring(G; degrees, basis=nothing,
                                  codimensions=nothing,
                                  basis_names=nothing,
                                  degree_convention=:cohomological)

Construct a presentation of the small rational quantum cohomology of `G` as
an Oscar graded quotient ring. The finite vector `degrees` specifies the
effective curve classes whose quantum corrections are included; degree zero
is included automatically.

The result is `(QH, classes, q)`. Here `QH` is the quotient, `classes[d + 1]`
contains the chosen codimension-`d` cohomology basis elements, and `q` contains
the Novikov generators for a basis of ``H_2(G;\mathbb Z)``.

For a Weyl-backed homogeneous graph, the polynomial Schubert basis is used
automatically. For the standard GKM graph of projective space, the basis is
``1,H,\ldots,H^n``. Other graphs must provide a polynomial GKM `basis`, its
`codimensions`, and optionally `basis_names`.

# Example
```jldoctest small_quantum_oscar
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta = GKMtools.curve_class(P1, Oscar.Edge(1, 2));

julia> QH, classes, q = small_quantum_cohomology_ring(P1; degrees=[beta]);

julia> h = only(classes[2]);

julia> h^2 == only(q)
true
```
"""
function small_quantum_cohomology_ring(
  G::AbstractGKMGraph;
  degrees::AbstractVector{<:CurveClass},
  basis=nothing,
  codimensions=nothing,
  basis_names=nothing,
  degree_convention::Symbol=:cohomological,
)
  degree_factor = degree_convention == :cohomological ? 2 :
    degree_convention == :codimension ? 1 : throw(ArgumentError(
      "degree_convention must be :cohomological or :codimension",
    ))

  polynomial_basis, basis_codimensions, names = _quantum_presentation_basis(
    G, basis, codimensions, basis_names,
  )
  zero_beta = _zero_curve_class(G)
  quantum_degrees = unique(vcat(CurveClass[zero_beta], collect(degrees)))
  foreach(beta -> _check_curve_class(G, beta), quantum_degrees)
  all(beta -> is_effective(G, beta), quantum_degrees) ||
    throw(ArgumentError("all quantum degrees must be effective"))

  H2 = GKM_second_homology(G).H2
  h2_rank = ngens(H2)
  q_names = ["q$i" for i in 1:h2_rank]
  q_degrees = [degree_factor * Int(chern_number(G, gen(H2, i))) for i in 1:h2_rank]
  any(<=(0), q_degrees) && throw(ArgumentError(
    "the chosen H2 generators must have positive first Chern degree",
  ))

  positive_positions = collect(2:length(polynomial_basis))
  variable_names = vcat(names[positive_positions], q_names)
  variable_degrees = vcat(
    degree_factor .* basis_codimensions[positive_positions], q_degrees,
  )
  uses_dummy_variable = isempty(variable_names)
  uses_dummy_variable && (variable_names = ["_zero"] ; variable_degrees = [degree_factor])
  raw_ring, _ = polynomial_ring(QQ, variable_names)
  graded_ring, variables = grade(raw_ring, variable_degrees)
  class_variables = variables[1:length(positive_positions)]
  q_variables = variables[length(positive_positions)+1:length(positive_positions)+h2_rank]

  relations = uses_dummy_variable ? [variables[1]] : elem_type(graded_ring)[]
  for i in positive_positions, j in i:length(polynomial_basis)
    rhs = zero(graded_ring)
    for beta in quantum_degrees
      product = quantum_product(G, beta, polynomial_basis[i], polynomial_basis[j])
      iszero(product) && continue
      coefficients = _ordinary_coefficients_in_basis(G, product, polynomial_basis)
      q_monomial = _novikov_monomial(q_variables, beta)
      for k in eachindex(coefficients)
        iszero(coefficients[k]) && continue
        basis_monomial = k == 1 ? one(graded_ring) : class_variables[k - 1]
        rhs += coefficients[k] * q_monomial * basis_monomial
      end
    end
    push!(relations, class_variables[i - 1] * class_variables[j - 1] - rhs)
  end

  quotient, quotient_map = quo(graded_ring, ideal(graded_ring, relations))
  flat_basis = vcat([one(quotient)], quotient_map.(class_variables))
  grouped_basis = [elem_type(quotient)[] for _ in 0:maximum(basis_codimensions)]
  for k in eachindex(flat_basis)
    push!(grouped_basis[basis_codimensions[k] + 1], flat_basis[k])
  end
  return quotient, grouped_basis, quotient_map.(q_variables)
end

function _quantum_presentation_basis(G, basis, codimensions, basis_names)
  if !isnothing(basis)
    classes = collect(basis)
    length(classes) == num_vertices(G) || throw(DimensionMismatch(
      "the cohomology basis must have $(num_vertices(G)) elements",
    ))
    all(c -> c isa GKMClass && graph(c) === G, classes) ||
      throw(ArgumentError("all basis elements must be GKM classes on G"))
    isnothing(codimensions) && throw(ArgumentError(
      "codimensions are required with an explicit basis",
    ))
    dims = Int.(codimensions)
    length(dims) == length(classes) || throw(DimensionMismatch(
      "codimensions and basis must have equal length",
    ))
    names = isnothing(basis_names) ? ["x$i" for i in eachindex(classes)] : String.(basis_names)
    length(names) == length(classes) || throw(DimensionMismatch(
      "basis_names and basis must have equal length",
    ))
    return classes, dims, names
  end
  return _automatic_quantum_presentation_basis(G)
end

function _automatic_quantum_presentation_basis(
  G::AbstractGKMGraph{R,V,F},
) where {R,V<:GeneralizedFlagVertex,F}
  order = sort(collect(1:num_vertices(G)); by=v -> (length(flag(vertices_structure(G)[v])), v))
  dims = [length(flag(vertices_structure(G)[v])) for v in order]
  classes = [schubert_basis(G, v; representation=:polynomial) for v in order]
  names = ["σ_" * label(G, v) for v in order]
  return classes, dims, names
end

function _automatic_quantum_presentation_basis(G::AbstractGKMGraph)
  n = num_vertices(G) - 1
  num_edges(G) == binomial(n + 1, 2) && rank_torus(G) == n + 1 ||
    throw(ArgumentError(
      "automatic quantum presentations are available for homogeneous Schubert graphs and standard projective spaces; provide basis and codimensions",
    ))
  t = collect(gens(equivariant_coefficient_ring(G)))
  H = polynomial_class(G, t)
  classes = [H^k for k in 0:n]
  return classes, collect(0:n), vcat(["1"], [k == 1 ? "h" : "h$k" for k in 1:n])
end

function _ordinary_coefficients_in_basis(G, product, basis)
  K = get_cohomology(G).localized_coefficient_ring
  B = matrix(K, length(basis), num_vertices(G), [K(c[v]) for c in basis for v in 1:num_vertices(G)])
  values = matrix(K, 1, num_vertices(G), K.(product.restrictions))
  coefficients = values * inv(B)
  origin = fill(QQ(0), ngens(equivariant_coefficient_ring(G)))
  return QQFieldElem[_specialize_quantum_coefficient(coefficients[k], origin) for k in 1:length(basis)]
end

function _specialize_quantum_coefficient(c, origin)
  numerator_value = evaluate(numerator(c), origin)
  denominator_value = evaluate(denominator(c), origin)
  iszero(denominator_value) && throw(ArgumentError(
    "the chosen equivariant basis does not specialize to an ordinary cohomology basis",
  ))
  return QQ(numerator_value) / QQ(denominator_value)
end

function _novikov_monomial(q_variables, beta::CurveClass)
  result = one(parent(first(q_variables)))
  for i in eachindex(q_variables)
    exponent = Int(beta[i])
    exponent >= 0 || throw(ArgumentError(
      "quantum degrees must have nonnegative coordinates in the H2 basis",
    ))
    result *= q_variables[i]^exponent
  end
  return result
end
