_zero_curve_class(G::AbstractGKMGraph) = zero(GKM_second_homology(G).H2)

function _check_curve_class(G::AbstractGKMGraph, beta::CurveClass)
  parent(beta) === GKM_second_homology(G).H2 ||
    throw(ArgumentError("the curve class belongs to a different second homology group"))
end

@doc raw"""
    small_equivariant_quantum_cohomology(G; degrees=[0])

Construct the small equivariant quantum cohomology of `G`, truncated to the
finite vector of effective curve classes `degrees`. The zero class is included
automatically.

# Example
```jldoctest small_quantum
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta = GKMtools.curve_class(P1, Oscar.Edge(1, 2));

julia> QH = GKMtools.small_equivariant_quantum_cohomology(P1; degrees=[beta]);

julia> length(truncation_degrees(QH))
2
```
"""
function small_equivariant_quantum_cohomology(G::AbstractGKMGraph;
    degrees::AbstractVector{<:CurveClass}=CurveClass[_zero_curve_class(G)])
  foreach(beta -> _check_curve_class(G, beta), degrees)
  all(beta -> is_effective(G, beta), degrees) ||
    throw(ArgumentError("all truncation degrees must be effective"))
  SmallEquivariantQuantumCohomology(
    G, unique(vcat(CurveClass[_zero_curve_class(G)], collect(degrees))))
end

function _coerce_quantum_coefficient(QH, value)
  G = graph(QH)
  if value isa GKMClass
    value.graph === G || throw(ArgumentError("coefficient belongs to a different graph"))
    return localize(value)
  elseif value isa AbstractVector
    return localized_class(G, value)
  end
  value * localize(unit_cohomology_ring(G))
end

@doc raw"""
    quantum_class(QH, value; degree=0)
    quantum_class(QH, coefficients)
    quantum_class(G, value; degree=0, degrees=[degree])

Construct a quantum class from a [`GKMClass`](@ref), a vector of fixed-point
restrictions, a scalar, or a dictionary of Novikov coefficients. Every degree
must belong to the truncation of `QH`.

# Example
```jldoctest small_quantum
julia> p = quantum_class(QH, point_class(P1, 1));

julia> q = quantum_class(QH, unit_cohomology_ring(P1); degree=beta);

julia> coefficient(p + q, beta) == localize(unit_cohomology_ring(P1))
true

julia> quantum_class(QH, [1, 1]) == one(QH)
true
```
"""
function quantum_class(QH::SmallEquivariantQuantumCohomology,
    values::AbstractDict{<:CurveClass})
  allowed = Set(QH.degrees)
  result = Dict{CurveClass,GKMClass}()
  for (beta, value) in values
    _check_curve_class(graph(QH), beta)
    beta in allowed || throw(ArgumentError("curve class $beta is outside the truncation"))
    c = _coerce_quantum_coefficient(QH, value)
    iszero(c) || (result[beta] = c)
  end
  SmallEquivariantQuantumClass(QH, result)
end

function quantum_class(QH::SmallEquivariantQuantumCohomology, value;
    degree::CurveClass=_zero_curve_class(graph(QH)))
  quantum_class(QH, Dict(degree => value))
end

function quantum_class(G::AbstractGKMGraph, value;
    degree::CurveClass=_zero_curve_class(G),
    degrees::AbstractVector{<:CurveClass}=CurveClass[degree])
  quantum_class(small_equivariant_quantum_cohomology(G; degrees), value; degree)
end

@doc raw"""
    QH_class(G, value; beta=nothing, degree=nothing, degrees=nothing)

Compatibility constructor for [`quantum_class`](@ref). The `beta` keyword is
an alias for `degree`.

# Example
```jldoctest qh_class
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta = GKMtools.curve_class(P1, Oscar.Edge(1, 2));

julia> q = QH_class(P1, 1; beta=beta);

julia> coefficient(q, beta) == localize(unit_cohomology_ring(P1))
true
```
"""
function QH_class(G::AbstractGKMGraph, value;
    beta::Union{Nothing,CurveClass}=nothing,
    degree::Union{Nothing,CurveClass}=nothing,
    degrees::Union{Nothing,AbstractVector{<:CurveClass}}=nothing)
  d = isnothing(degree) ? (isnothing(beta) ? _zero_curve_class(G) : beta) : degree
  ds = isnothing(degrees) ? CurveClass[d] : degrees
  quantum_class(G, value; degree=d, degrees=ds)
end

@doc raw"""
    quantum_coefficients(c)

Return a copy of the dictionary of nonzero Novikov coefficients of `c`.
"""
quantum_coefficients(c::SmallEquivariantQuantumClass) = copy(c.coefficients)

@doc raw"""
    coefficient(c, beta)

Return the coefficient of ``q^\beta`` in `c`. A degree outside the support has
the zero GKM class as coefficient.

# Example
```jldoctest quantum_coefficients
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta = GKMtools.curve_class(P1, Oscar.Edge(1, 2));

julia> QH = GKMtools.small_equivariant_quantum_cohomology(P1; degrees=[beta]);

julia> c = quantum_class(QH, Dict(beta => [1, 1]));

julia> length(quantum_coefficients(c))
1

julia> iszero(coefficient(c, zero(GKM_second_homology(P1).H2)))
true
```
"""
function coefficient(c::SmallEquivariantQuantumClass, beta::CurveClass)
  _check_curve_class(graph(c), beta)
  get(c.coefficients, beta, zero(localize(unit_cohomology_ring(graph(c)))))
end

function _check_same_quantum_parent(a, b)
  parent(a) === parent(b) || throw(ArgumentError("quantum classes have different parents"))
end

Base.zero(QH::SmallEquivariantQuantumCohomology) =
  SmallEquivariantQuantumClass(QH, Dict{CurveClass,GKMClass}())
Base.zero(c::SmallEquivariantQuantumClass) = zero(parent(c))
Base.one(QH::SmallEquivariantQuantumCohomology) =
  quantum_class(QH, unit_cohomology_ring(graph(QH)))
Base.one(c::SmallEquivariantQuantumClass) = one(parent(c))
Base.iszero(c::SmallEquivariantQuantumClass) = isempty(c.coefficients)
Base.isone(c::SmallEquivariantQuantumClass) = c == one(c)

function Base.:+(a::SmallEquivariantQuantumClass, b::SmallEquivariantQuantumClass)
  _check_same_quantum_parent(a, b)
  values = copy(a.coefficients)
  for (beta, value) in b.coefficients
    values[beta] = get(values, beta, zero(value)) + value
    iszero(values[beta]) && delete!(values, beta)
  end
  SmallEquivariantQuantumClass(parent(a), values)
end
Base.:-(c::SmallEquivariantQuantumClass) = (-1) * c
Base.:-(a::SmallEquivariantQuantumClass, b::SmallEquivariantQuantumClass) = a + (-b)

function Base.:*(scalar, c::SmallEquivariantQuantumClass)
  values = Dict{CurveClass,GKMClass}()
  for (beta, value) in c.coefficients
    product = scalar * value
    iszero(product) || (values[beta] = product)
  end
  SmallEquivariantQuantumClass(parent(c), values)
end
Base.:*(c::SmallEquivariantQuantumClass, scalar) = scalar * c

function _homogeneous_degree(c::GKMClass)
  degrees = Int[]
  for value in c.restrictions
    iszero(value) && continue
    _is_homogeneous(value) || return nothing
    push!(degrees, _get_degree(value))
  end
  isempty(degrees) && return nothing
  all(==(first(degrees)), degrees) || return nothing
  return first(degrees)
end

function _quantum_product_degree(G, beta, a, b)
  degree_a = _homogeneous_degree(a)
  degree_b = _homogeneous_degree(b)
  (isnothing(degree_a) || isnothing(degree_b)) && return nothing
  return degree_a + degree_b - Int(chern_number(G, beta))
end

_quantum_product_vanishes_by_degree(G, beta, a, b) =
  something(_quantum_product_degree(G, beta, a, b), 0) < 0

@doc raw"""
    quantum_product(G, beta, a, b; show_bar=false)

Return the coefficient of `q^beta` in the small equivariant quantum product
of the GKM classes `a` and `b`. Degree zero is the ordinary cup product.
Positive effective degrees are computed using [`gromov_witten_nomarks`](@ref),
with the two factors and each fixed-point class integrated over the image
curve.

# Example
```jldoctest quantum_product
julia> P1 = GKMtools.projective_space(GKMtools.GKMGraph, 1);

julia> beta0 = zero(GKM_second_homology(P1).H2);

julia> p = point_class(P1, 1);

julia> quantum_product(P1, beta0, p, p) == localize(p^2)
true
```
"""
function quantum_product(G::AbstractGKMGraph, beta::CurveClass,
    a::GKMClass, b::GKMClass; show_bar::Bool=false)
  _check_curve_class(G, beta)
  a.graph === G && b.graph === G ||
    throw(ArgumentError("the factors must belong to the given GKM graph"))
  beta == _zero_curve_class(G) && return localize(a * b)
  is_effective(G, beta) ||
    return zero(localize(unit_cohomology_ring(G)))
  _quantum_product_vanishes_by_degree(G, beta, a, b) &&
    return zero(localize(unit_cohomology_ring(G)))
  if _quantum_product_degree(G, beta, a, b) == 0
    invariant = gromov_witten_nomarks(
      G, beta, GKMClass[a, b, point_class(G, 1)]; show_bar,
    )
    return invariant * localize(unit_cohomology_ring(G))
  end
  class_products = [
    GKMClass[a, b, point_class(G, v)]
    for v in 1:num_vertices(G)
  ]
  localized_class(G, gromov_witten_nomarks(G, beta, class_products; show_bar))
end

function Base.:*(a::SmallEquivariantQuantumClass, b::SmallEquivariantQuantumClass)
  _check_same_quantum_parent(a, b)
  QH = parent(a)
  values = Dict{CurveClass,GKMClass}()
  for gamma in QH.degrees
    value = zero(localize(unit_cohomology_ring(graph(QH))))
    for (beta1, a1) in a.coefficients, (beta2, b2) in b.coefficients
      delta = gamma - beta1 - beta2
      is_effective(graph(QH), delta) || continue
      value += quantum_product(graph(QH), delta, a1, b2)
    end
    iszero(value) || (values[gamma] = value)
  end
  SmallEquivariantQuantumClass(QH, values)
end

function Base.:^(c::SmallEquivariantQuantumClass, n::Integer)
  n >= 0 || throw(DomainError(n, "a quantum class cannot have a negative power"))
  result, factor, exponent = one(c), c, n
  while exponent > 0
    isodd(exponent) && (result *= factor)
    exponent >>= 1
    exponent > 0 && (factor *= factor)
  end
  result
end

Base.:(==)(a::SmallEquivariantQuantumClass, b::SmallEquivariantQuantumClass) =
  parent(a) === parent(b) && a.coefficients == b.coefficients

Base.:+(c::SmallEquivariantQuantumClass, value) = c + quantum_class(parent(c), value)
Base.:+(value, c::SmallEquivariantQuantumClass) = c + value
Base.:-(c::SmallEquivariantQuantumClass, value) = c + (-value)
Base.:-(value, c::SmallEquivariantQuantumClass) = quantum_class(parent(c), value) - c
