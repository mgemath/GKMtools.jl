###############################################################################
#
# Wall--Jupp invariants of signed 3-valent GKM graphs
#
###############################################################################

"""
An integral Wall--Jupp--Zubr invariant system obtained from a signed 3-valent GKM
graph. Coordinates in `H2` are columns. `p1` and `c2` are row covectors, and
`mu_packed` stores one value for every triple in `mu_triples`.

The construction records, but cannot verify from combinatorial data alone, the
geometric assumptions that the graph comes from a closed simply connected
six-manifold with torsion-free integral cohomology.
"""
struct _GKMInvariantSystem
  H2_rank::Int
  H2::Any
  basis_localizations::ZZMatrix
  mu_packed::ZZMatrix
  mu_triples::Vector{NTuple{3, Int}}
  cubic_ring::Any
  cubic::Any
  c1::ZZMatrix
  w2::Any
  p1::ZZMatrix
  c2::ZZMatrix
  b3::ZZRingElem
  euler::ZZRingElem
  characteristic_numbers::NamedTuple
  root::Int
  weight_rank::Int
  diagnostics::Dict{Symbol, Any}
  cache::Dict{Symbol, Any}
end

"""Result of `compare_systems`; positive answers carry a verified matrix."""
struct _SystemComparisonResult
  status::Symbol
  witness::Any
  obstruction::Any
  diagnostics::Vector{Any}
  function _SystemComparisonResult(status, witness, obstruction, diagnostics)
    status in (:equivalent, :not_equivalent, :unknown) ||
      throw(ArgumentError("invalid comparison status $status"))
    status == :equivalent && !(witness isa ZZMatrix) &&
      throw(ArgumentError("an :equivalent result requires an integral witness"))
    status == :not_equivalent && isnothing(obstruction) &&
      throw(ArgumentError("a :not_equivalent result requires an obstruction certificate"))
    status == :unknown && (!isnothing(witness) || !isnothing(obstruction)) &&
      throw(ArgumentError("an :unknown result cannot carry a witness or obstruction"))
    return new(status, witness, obstruction, diagnostics)
  end
end

function Base.show(io::IO, result::_SystemComparisonResult)
  print(io, "System comparison: ", result.status)
  !isnothing(result.obstruction) && print(io, " (", result.obstruction, ")")
end

function Base.show(io::IO, S::_GKMInvariantSystem)
  print(
    io,
    "6-dimensional invariant system (rank(H^2) = ", S.H2_rank,
    ", b3 = ", S.b3,
    ", Euler characteristic = ", S.euler,
    ")",
  )
end

_root_index(G::AbstractGKM_graph, root::Integer) = Int(root)
function _root_index(G::AbstractGKM_graph, root::AbstractString)
  i = findfirst(==(root), G.labels)
  @req !isnothing(i) "Unknown root vertex label $root."
  return i
end

function _weight_column(G::AbstractGKM_graph, e::Edge)
  a = _w(G, e)
  return matrix(ZZ, rank_torus(G), 1, [ZZ(a[i]) for i in 1:rank_torus(G)])
end

function _weight_column(G::AbstractGKM_graph, v::Int, flag::Int)
  a = G.weights_at_vertex[v][flag]
  return matrix(ZZ, rank_torus(G), 1, [ZZ(a[i]) for i in 1:rank_torus(G)])
end

function _integral_multiple(a::ZZMatrix, b::ZZMatrix)
  # Decide whether a=q*b for q in ZZ. The caller guarantees b is nonzero.
  j = findfirst(i -> !iszero(b[i, 1]), 1:nrows(b))
  isnothing(j) && return false
  q = QQ(a[j, 1]) // QQ(b[j, 1])
  denominator(q) == 1 || return false
  return all(i -> QQ(a[i, 1]) == q * QQ(b[i, 1]), 1:nrows(a))
end

"""
    _validate_signed_gkm_graph(G; check_characteristic_classes=true)

Validate the signed, integral, compact 3-valent input. This begins with
GKMtools' `isvalid(G)` and adds the dimension-six and integral characteristic
class checks needed by `system_of_invariants_6d`. The returned named tuple has
fields `valid`, `errors`, `warnings`, `effective_rank`, and `rank_defect`.
"""
function _validate_signed_gkm_graph(
  G::AbstractGKM_graph;
  check_characteristic_classes::Bool=true,
)
  errors = String[]
  warnings = String[]

  isvalid(G; printDiagnostics=false) || push!(errors, "GKMtools.isvalid reports an invalid graph")
  base_ring(G.M) == ZZ || push!(errors, "the weight lattice must be a free ZZ-module")
  is_compact(G) || push!(errors, "the graph must be compact (no standalone flags)")
  is_connected(G.g) || push!(errors, "the underlying graph must be connected")
  edge_pairs = [(min(src(e), dst(e)), max(src(e), dst(e))) for e in edges(G.g)]
  all(e -> src(e) != dst(e), edges(G.g)) || push!(errors, "the underlying graph must have no loops")
  length(unique(edge_pairs)) == length(edge_pairs) || push!(errors, "the underlying graph must be simple")
  valency(G) == 3 || push!(errors, "the graph must be 3-valent")

  k = rank_torus(G)
  all_weights = ZZMatrix[]
  for v in 1:n_vertices(G.g), i in 1:length(G.weights_at_vertex[v])
    a = try
      _weight_column(G, v, i)
    catch
      push!(errors, "weight $i at vertex $v is not integral")
      continue
    end
    iszero(a) && push!(errors, "weight $i at vertex $v is zero")
    push!(all_weights, a)
  end

  effective_rank = 0
  if !isempty(all_weights)
    W = zero_matrix(ZZ, k, length(all_weights))
    for (j, a) in enumerate(all_weights), i in 1:k
      W[i, j] = a[i, 1]
    end
    effective_rank = rank(change_base_ring(QQ, W))
  end
  rank_defect = k - effective_rank
  rank_defect > 0 && push!(warnings, "the input torus has rational ineffective rank $rank_defect")

  nonprimitive_edges = Edge[]
  for e in edges(G.g)
    alpha = _weight_column(G, e)
    gcd_entries = foldl(gcd, (abs(alpha[i, 1]) for i in 1:k); init=ZZ(0))
    gcd_entries > 1 && push!(nonprimitive_edges, e)
  end
  if !isempty(nonprimitive_edges)
    push!(
      warnings,
      "$(length(nonprimitive_edges)) edge weight(s) are non-primitive; identifying integral graph cohomology with integral equivariant cohomology may require additional isotropy or Chang--Skjelbred hypotheses",
    )
  end

  if check_characteristic_classes && isempty(errors)
    # Also use the package's native GKM-class check.
    c1T = first_chern_class(G)
    p1T = first_pontryagin_class(G)
    is_gkm_class(c1T, G) || push!(errors, "the tangent c1 assignment is not a GKM class")
    is_gkm_class(p1T, G) || push!(errors, "the tangent p1 assignment is not a GKM class")

    Rz, z = polynomial_ring(ZZ, ["u$i" for i in 1:k])
    linear(a) = sum(a[i, 1] * z[i] for i in 1:k; init=zero(Rz))
    c1_vectors = ZZMatrix[]
    p1_polys = Any[]
    for v in 1:n_vertices(G.g)
      av = [_weight_column(G, v, i) for i in 1:3]
      push!(c1_vectors, sum(av; init=zero_matrix(ZZ, k, 1)))
      push!(p1_polys, sum(linear(a)^2 for a in av; init=zero(Rz)))
    end
    for e in edges(G.g)
      u, v = src(e), dst(e)
      alpha = _weight_column(G, e)
      _integral_multiple(c1_vectors[u] - c1_vectors[v], alpha) ||
        push!(errors, "c1 fails integral edge divisibility along $e")
      divides(p1_polys[u] - p1_polys[v], linear(alpha))[1] ||
        push!(errors, "p1 fails integral polynomial divisibility along $e")
    end
  end

  return (
    valid=isempty(errors), errors=errors, warnings=warnings,
    effective_rank=effective_rank, rank_defect=rank_defect,
  )
end

"""
    _build_degree_two_relation_matrix(G, root=1)

Construct the integral relation matrix for normalized degree-two GKM classes.
The unknowns are the `k` coordinates at every non-root vertex followed by one
integer edge quotient for each edge.
"""
function _build_degree_two_relation_matrix(
  G::AbstractGKM_graph,
  root=first(vertices(G.g)),
)
  rho = _root_index(G, root)
  n = n_vertices(G.g)
  @req 1 <= rho <= n "Root vertex is out of range."
  k = rank_torus(G)
  edge_list = collect(edges(G.g))
  vertex_offset = Dict(v => k * (j - 1) for (j, v) in enumerate(filter(!=(rho), 1:n)))
  nf = k * (n - 1)
  M = zero_matrix(ZZ, k * length(edge_list), nf + length(edge_list))

  for (j, e) in enumerate(edge_list)
    u, v = src(e), dst(e)
    alpha = _weight_column(G, e)
    for i in 1:k
      row = k * (j - 1) + i
      u != rho && (M[row, vertex_offset[u] + i] = -1)
      v != rho && (M[row, vertex_offset[v] + i] = 1)
      M[row, nf + j] = -alpha[i, 1]
    end
  end
  return M
end

"""
    _integral_kernel_basis(M)

Return a saturated integral basis, in columns, of the right kernel of the
integer matrix `M`. Smith transforms are used so saturation does not depend on
clearing denominators in a rational nullspace.
"""
function _integral_kernel_basis(M::ZZMatrix)
  D, U, V = snf_with_transform(M)
  @assert U * M * V == D
  rk = rank(change_base_ring(QQ, M))
  K = rk == ncols(M) ? zero_matrix(ZZ, ncols(M), 0) : V[:, (rk + 1):ncols(M)]
  @assert M * K == zero_matrix(ZZ, nrows(M), ncols(K))
  @assert rank(change_base_ring(QQ, K)) == ncols(M) - rk
  return K
end

"""
    _normalized_h2_basis(G, root=1)

Return the `n*k` by `b2` localization matrix whose columns form the saturated
integral lattice of ordinary degree-two classes in the root-zero gauge.
"""
function _normalized_h2_basis(G::AbstractGKM_graph, root=first(vertices(G.g)))
  rho = _root_index(G, root)
  n, k = n_vertices(G.g), rank_torus(G)
  M = _build_degree_two_relation_matrix(G, rho)
  K = _integral_kernel_basis(M)
  nf = k * (n - 1)
  P = K[1:nf, :]
  L = zero_matrix(ZZ, n * k, ncols(P))
  source_block = 0
  for v in 1:n
    v == rho && continue
    for i in 1:k, j in 1:ncols(P)
      L[k * (v - 1) + i, j] = P[k * source_block + i, j]
    end
    source_block += 1
  end
  @assert rank(change_base_ring(QQ, L)) == ncols(L)
  return L
end

"""Deterministic packed ordering for a symmetric trilinear tensor of rank `r`."""
_symmetric_triples(r::Integer) = [(i, j, k) for i in 1:r for j in i:r for k in j:r]

function _triple_index(S::_GKMInvariantSystem, i::Int, j::Int, k::Int)
  t = Tuple(sort([i, j, k]))
  index = get!(S.cache, :triple_index) do
    Dict(triple => position for (position, triple) in enumerate(S.mu_triples))
  end
  p = get(index, t, nothing)
  @req !isnothing(p) "Tensor index is out of range."
  return p
end

"""Return `mu(e_i,e_j,e_k)` from the full packed tensor (indices are sorted)."""
_mu_entry(S::_GKMInvariantSystem, i::Int, j::Int, k::Int) =
  S.mu_packed[1, _triple_index(S, i, j, k)]

function _vector_entry(x, i)
  x isa AbstractMatrix && return x[i, 1]
  return x[i]
end

"""Evaluate the full symmetric trilinear tensor on three column vectors."""
function _mu_eval(S::_GKMInvariantSystem, x, y, z)
  r = S.H2_rank
  @req length(x) == r && length(y) == r && length(z) == r "Vector has the wrong rank."
  r == 0 && return ZZ(0)
  sample = S.mu_packed[1, 1] * _vector_entry(x, 1) * _vector_entry(y, 1) * _vector_entry(z, 1)
  result = zero(sample)
  for (q, (i, j, k)) in enumerate(S.mu_triples)
    mu = S.mu_packed[1, q]
    xi, xj, xk = _vector_entry(x, i), _vector_entry(x, j), _vector_entry(x, k)
    yi, yj, yk = _vector_entry(y, i), _vector_entry(y, j), _vector_entry(y, k)
    zi, zj, zk = _vector_entry(z, i), _vector_entry(z, j), _vector_entry(z, k)
    if i == k
      result += mu * xi * yi * zi
    elseif i == j
      result += mu * (xi * yi * zk + xi * yk * zi + xk * yi * zi)
    elseif j == k
      result += mu * (xi * yj * zj + xj * yi * zj + xj * yj * zi)
    else
      result += mu * (
        xi * yj * zk + xi * yk * zj +
        xj * yi * zk + xj * yk * zi +
        xk * yi * zj + xk * yj * zi
      )
    end
  end
  return result
end

"""Return the matrix of `Sym^2(H) -> H^*` defined by the trilinear tensor."""
function _mu_flattening(S::_GKMInvariantSystem; base_ring=ZZ)
  r = S.H2_rank
  pairs = [(j, k) for j in 1:r for k in j:r]
  F = zero_matrix(base_ring, r, length(pairs))
  for i in 1:r, (q, (j, k)) in enumerate(pairs)
    F[i, q] = base_ring(_mu_entry(S, i, j, k))
  end
  return F
end

"""Return `B_x`, where `B_x[i,j]=mu(x,e_i,e_j)`."""
function _contraction_matrix(S::_GKMInvariantSystem, x; base_ring=ZZ)
  r = S.H2_rank
  B = zero_matrix(base_ring, r, r)
  for i in 1:r, j in i:r
    value = sum(base_ring(_vector_entry(x, a)) * base_ring(_mu_entry(S, a, i, j)) for a in 1:r; init=zero(base_ring))
    B[i, j] = B[j, i] = value
  end
  return B
end

function _multiplicity(i, j, k) # This assumes i <= j <= k.
  i == k && return 1
  (i == j || j == k) && return 3
  return 6
end

"""Reconstruct `mu(z,z,z)` over a requested coefficient ring."""
function _diagonal_cubic(S::_GKMInvariantSystem; base_ring=ZZ)
  R, z = polynomial_ring(base_ring, ["z$i" for i in 1:S.H2_rank])
  F = zero(R)
  for (q, (i, j, k)) in enumerate(S.mu_triples)
    F += _multiplicity(i, j, k) * base_ring(S.mu_packed[1, q]) * z[i] * z[j] * z[k]
  end
  return R, F
end

"""Return `det(sum_a z_a B_ea)`, the determinant polynomial of contractions."""
function _determinant_polynomial(S::_GKMInvariantSystem; base_ring=QQ)
  R, z = polynomial_ring(base_ring, ["z$i" for i in 1:S.H2_rank])
  r = S.H2_rank
  B = zero_matrix(R, r, r)
  for i in 1:r, j in 1:r
    B[i, j] = sum(base_ring(_mu_entry(S, a, i, j)) * z[a] for a in 1:r; init=zero(R))
  end
  return R, det(B)
end

function _factorization_signature(f)
  iszero(f) && return (zero=true, total_degree=-1, factors=Tuple{Int, Int}[])
  factors = Tuple{Int, Int}[]
  for (factor_polynomial, multiplicity) in factor(f)
    push!(factors, (total_degree(factor_polynomial), Int(multiplicity)))
  end
  sort!(factors)
  return (zero=false, total_degree=total_degree(f), factors=factors)
end

function _rational_polynomial_signatures(S::_GKMInvariantSystem)
  return get!(S.cache, :rational_polynomial_signatures) do
    _, cubic = _diagonal_cubic(S; base_ring=QQ)
    _, determinant = _determinant_polynomial(S; base_ring=QQ)
    (
      cubic=_factorization_signature(cubic),
      contraction_determinant=_factorization_signature(determinant),
    )
  end
end

function _rational_polynomial_obstruction(S1::_GKMInvariantSystem, S2::_GKMInvariantSystem)
  signatures1 = _rational_polynomial_signatures(S1)
  signatures2 = _rational_polynomial_signatures(S2)
  signatures1.cubic == signatures2.cubic ||
    return "diagonal cubics have different rational factorization data"
  signatures1.contraction_determinant == signatures2.contraction_determinant ||
    return "contraction determinant polynomials have different rational factorization data"
  return nothing
end

function _class_from_localization(G::AbstractGKM_graph, L::ZZMatrix, column::Int)
  R = G.equivariantCohomology
  t = gens(R.coeffRing)
  k = rank_torus(G)
  values = QQMPolyRingElem[]
  for v in 1:n_vertices(G.g)
    p = zero(R.coeffRing)
    for i in 1:k
      p += L[k * (v - 1) + i, column] * t[i]
    end
    push!(values, p)
  end
  return R.cohomRing(values)
end

function _extract_integral_constant(q; context="localized integral")
  numerator_q, denominator_q = numerator(q), denominator(q)
  @req is_constant(numerator_q) && is_constant(denominator_q) "$context is not constant: $q"
  exponent = zeros(Int, length(gens(parent(numerator_q))))
  value = coeff(numerator_q, exponent) // coeff(denominator_q, exponent)
  @req denominator(value) == 1 "$context is not integral: $value"
  return ZZ(numerator(value))
end

function _linear_coefficients(p, k::Int)
  @req total_degree(p) <= 1 "Expected a linear localization, got $p."
  result = zero_matrix(ZZ, k, 1)
  for i in 1:k
    exponent = zeros(Int, k)
    exponent[i] = 1
    c = coeff(p, exponent)
    @req denominator(c) == 1 "Linear localization has a nonintegral coefficient."
    result[i, 1] = ZZ(numerator(c))
  end
  return result
end

@doc raw"""
    system_of_invariants_6d(G; root=first(vertices(G.g)), check=true)

Construct the Wall--Jupp--Žubr [MR215313, MR314074, MR970082](@cite) _system of invariants_ of the compact Hamiltonian GKM space
$X$ of real dimension six encoded by the GKM graph $G$.
The classification interpretation assumes that $X$ is closed, simply connected,
and oriented, that its integral cohomology is torsion-free, and that the signed
graph computes its integral equivariant cohomology.

The result has fields for

- the free $\mathbb{Z}$-module $H^2(X;\mathbb{Z})$,
- the odd Betti number $b_3$,
- the [`first_pontryagin_class`](@ref) $p_1(T_X)$,
- the second Stiefel--Whitney class $w_2(T_X)$,
- the tensor $\mu_X\colon H\times H \times H \rightarrow \mathbb{Z}$ given by $\mu_X(a,b,c) = \int_X abc$.

It also includes $c_1(T_X)$ which encodes the homotopy class of the almost complex structure.
In the isolated-fixed-point GKM setting used here, odd cohomology vanishes, so
$b_3=0$; systems for connected sums with $S^3\times S^3$ are not represented.
The underlying `GKM_graph` type is simple, so parallel-edge GKM multigraphs are
not supported. Non-primitive weights are retained verbatim and recorded in the
validation warnings.

# Optional arguments:
- `root`: the vertex at which degree-two localizations are normalized to zero.
- `check::Bool`: `true` by default. Run all graph and characteristic-class
  validation. With `false`, these potentially expensive validation passes are
  skipped; the algebraic preconditions needed by the construction are still
  enforced.
"""
function system_of_invariants_6d(
  G::AbstractGKM_graph;
  root=first(vertices(G.g)),
  check::Bool=true,
)
  rho = _root_index(G, root)
  validation = if check
    checked = _validate_signed_gkm_graph(G)
    @req checked.valid join(checked.errors, "; ")
    checked
  else
    (
      valid=true,
      errors=String[],
      warnings=["graph validation was skipped because check=false"],
      effective_rank=nothing,
      rank_defect=nothing,
    )
  end
  @req base_ring(G.M) == ZZ "Integral invariants require a ZZ weight lattice."
  @req valency(G) == 3 "The input must be 3-valent."

  n, k = n_vertices(G.g), rank_torus(G)
  relation_matrix = _build_degree_two_relation_matrix(G, rho)
  rational_nullity = ncols(relation_matrix) - rank(change_base_ring(QQ, relation_matrix))
  L = _normalized_h2_basis(G, rho)
  r = ncols(L)
  @assert r == rational_nullity

  betti = betti_numbers(G)
  @req length(betti) == 4 "A 3-valent graph must have four even Betti entries."
  @req r == betti[2] "The integral degree-two lattice rank $r disagrees with betti_numbers(G)=$(betti[2])."
  b3_int = 2 + 2 * r - n
  @req b3_int >= 0 && iseven(b3_int) "Poincare duality gives an invalid b3=$b3_int."
  @req b3_int == 0 "Isolated-fixed-point GKM graphs must have b3=0; connected-sum S^3 x S^3 contributions are unsupported."

  h = [_class_from_localization(G, L, i) for i in 1:r]
  check && @assert all(is_gkm_class(x, G) for x in h)

  triples = _symmetric_triples(r)
  MU = zero_matrix(ZZ, 1, length(triples))
  for (q, (i, j, ell)) in enumerate(triples)
    MU[1, q] = _extract_integral_constant(integrate(h[i] * h[j] * h[ell], G); context="mu[$i,$j,$ell]")
  end

  c1T = first_chern_class(G)
  C1loc = zero_matrix(ZZ, n * k, 1)
  for v in 1:n
    block = _linear_coefficients(c1T[v] - c1T[rho], k)
    for i in 1:k
      C1loc[k * (v - 1) + i, 1] = block[i, 1]
    end
  end
  c1_solvable, c1 = can_solve_with_solution(L, C1loc; side=:right)
  @req c1_solvable "The normalized first Chern class is not in the saturated H2 lattice."
  @assert L * c1 == C1loc
  w2 = change_base_ring(GF(2), c1)

  p1T = first_pontryagin_class(G)
  p1 = zero_matrix(ZZ, 1, r)
  c2_localized = zero_matrix(ZZ, 1, r)
  c2T = chern_class(G, 2)
  for i in 1:r
    p1[1, i] = _extract_integral_constant(integrate(p1T * h[i], G); context="p1[$i]")
    c2_localized[1, i] = _extract_integral_constant(integrate(c2T * h[i], G); context="c2[$i]")
  end

  # Create a temporary tensor-only view for coordinate contractions.
  Rz, z = polynomial_ring(ZZ, ["z$i" for i in 1:r])
  cubic = zero(Rz)
  for (q, (i, j, ell)) in enumerate(triples)
    cubic += _multiplicity(i, j, ell) * MU[1, q] * z[i] * z[j] * z[ell]
  end
  # Use an explicit full sum: diagonal-polynomial polarization is invalid in
  # characteristics 2 and 3 and is deliberately not used anywhere here.
  full_mu(x, y, zvec) = sum(
    MU[1, _triple_index_from_list(triples, i, j, ell)] * x[i, 1] * y[j, 1] * zvec[ell, 1]
    for i in 1:r for j in 1:r for ell in 1:r; init=ZZ(0)
  )

  l11 = zero_matrix(ZZ, 1, r)
  for i in 1:r
    ei = zero_matrix(ZZ, r, 1)
    ei[i, 1] = 1
    l11[1, i] = full_mu(c1, c1, ei)
  end
  c2 = zero_matrix(ZZ, 1, r)
  for i in 1:r
    d = l11[1, i] - p1[1, i]
    @req iseven(d) "(c1^2-p1)/2 is nonintegral on basis vector $i."
    c2[1, i] = divexact(d, ZZ(2))
  end
  @req c2 == c2_localized "The localized c2 disagrees with (c1^2-p1)/2."

  c1_cubed = full_mu(c1, c1, c1)
  p1_c1 = (p1 * c1)[1, 1]
  c1_c2 = (c2 * c1)[1, 1]
  @assert 2 * c1_c2 == c1_cubed - p1_c1
  c3 = _extract_integral_constant(integrate(chern_class(G, 3), G); context="c3")
  @req c3 == n "Localized c3=$c3 differs from the fixed-point Euler characteristic $n."

  diagnostics = Dict{Symbol, Any}(
    :validation => validation,
    :assumptions => [
      "the graph is realized by a closed simply connected oriented six-manifold",
      "integral equivariant graph cohomology computes its torsion-free integral cohomology",
      "the signed weights are its actual complex tangent weights",
    ],
    :rational_nullity => rational_nullity,
    :betti_numbers => betti,
  )
  characteristic = (c1_cubed=c1_cubed, c1_c2=c1_c2, c3=ZZ(n), p1_c1=p1_c1)
  S = _GKMInvariantSystem(
    r, free_module(ZZ, r), L, MU, triples, Rz, cubic, c1, w2, p1, c2,
    ZZ(b3_int), ZZ(n), characteristic, rho, k, diagnostics, Dict{Symbol, Any}(),
  )
  realizability = _realizability_diagnostics(S)
  diagnostics[:realizability] = realizability
  @req realizability.realizable "The computed system fails the Wall--Jupp Wu or Riemann--Roch realizability congruences."
  @assert _validate_invariant_system(S)
  return S
end

function _triple_index_from_list(triples, i, j, k)
  t = Tuple(sort([i, j, k]))
  p = findfirst(==(t), triples)
  @assert !isnothing(p)
  return p
end


"""Return the stored characteristic-number named tuple."""
_characteristic_numbers(S::_GKMInvariantSystem) = S.characteristic_numbers

function _multiindices_total_degree(r::Int, degree::Int)
  result = Vector{Vector{Int}}()
  current = zeros(Int, r)
  function extend(position, remaining)
    if position > r
      push!(result, copy(current))
      return
    end
    for value in 0:remaining
      current[position] = value
      extend(position + 1, remaining - value)
    end
  end
  extend(1, degree)
  return result
end

function _riemann_roch_index_numerator(S::_GKMInvariantSystem, x)
  return 2 * _mu_eval(S, x, x, x) +
    3 * _mu_eval(S, S.c1, x, x) +
    _mu_eval(S, S.c1, S.c1, x) +
    sum(S.c2[1, i] * _vector_entry(x, i) for i in 1:S.H2_rank; init=ZZ(0))
end

function _realizability_diagnostics(S::_GKMInvariantSystem)
  r = S.H2_rank
  wu_failures = Tuple{Int, Int}[]
  for i in 1:r, j in 1:r
    rhs = sum(
      Int(lift(ZZ, S.w2[a, 1])) * Int(mod(_mu_entry(S, a, i, j), 2))
      for a in 1:r; init=0,
    )
    lhs = Int(mod(_mu_entry(S, i, i, j) + _mu_entry(S, i, j, j), 2))
    mod(lhs - rhs, 2) == 0 || push!(wu_failures, (i, j))
  end

  rr_failures = Vector{Vector{Int}}()
  # The coefficients in the multivariate binomial basis are the iterated
  # forward differences at zero. Divisibility of all coefficients by 12 is
  # equivalent to the Riemann--Roch numerator being divisible by 12 for every
  # integral class, and requires only O(r^3) coefficient checks.
  for alpha in _multiindices_total_degree(r, 3)
    coefficient = ZZ(0)
    ranges = ntuple(i -> 0:alpha[i], r)
    for beta_tuple in Iterators.product(ranges...)
      beta = collect(beta_tuple)
      weight = prod(binomial(alpha[i], beta[i]) for i in 1:r; init=1)
      isodd(sum(alpha) - sum(beta)) && (weight = -weight)
      coefficient += weight * _riemann_roch_index_numerator(S, beta)
    end
    mod(coefficient, 12) == 0 || push!(rr_failures, alpha)
  end
  return (
    realizable=isempty(wu_failures) && isempty(rr_failures),
    wu_failures=wu_failures,
    riemann_roch_failures=rr_failures,
  )
end

"""Check internal dimensions and characteristic-class identities exactly."""
function _validate_invariant_system(S::_GKMInvariantSystem)
  r = S.H2_rank
  size(S.c1) == (r, 1) || return false
  size(S.w2) == (r, 1) || return false
  size(S.p1) == (1, r) || return false
  size(S.c2) == (1, r) || return false
  size(S.mu_packed) == (1, binomial(r + 2, 3)) || return false
  S.mu_triples == _symmetric_triples(r) || return false
  S.w2 == change_base_ring(GF(2), S.c1) || return false
  for i in 1:r
    ei = zero_matrix(ZZ, r, 1)
    ei[i, 1] = 1
    _mu_eval(S, S.c1, S.c1, ei) - S.p1[1, i] == 2 * S.c2[1, i] || return false
  end
  _mu_eval(S, S.c1, S.c1, S.c1) == S.characteristic_numbers.c1_cubed || return false
  (S.c2 * S.c1)[1, 1] == S.characteristic_numbers.c1_c2 || return false
  S.euler == S.characteristic_numbers.c3 || return false
  _realizability_diagnostics(S).realizable || return false
  return true
end

function _require_valid_invariant_system(S::_GKMInvariantSystem, label::AbstractString)
  valid = try
    _validate_invariant_system(S)
  catch error
    throw(ArgumentError(
      "$label input is malformed: $(sprint(showerror, error))",
    ))
  end
  valid || throw(ArgumentError(
    "$label input is malformed or fails the Wall--Jupp realizability congruences",
  ))
  return nothing
end

"""
    _verify_system_isomorphism(S1, S2, A;
        preserve_almost_complex=false)

Verify exactly that the integral matrix `A` maps column coordinates of `S1` to
column coordinates of `S2`. A determinant of either sign is allowed on `H^2`.
When `preserve_almost_complex=true`, `A` must additionally preserve `c1`.
"""
function _verify_system_isomorphism(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem,
  A::ZZMatrix;
  preserve_almost_complex::Bool=false,
)
  r = S1.H2_rank
  S2.H2_rank == r || return false
  size(A) == (r, r) || return false
  abs(det(A)) == 1 || return false

  preserve_almost_complex && A * S1.c1 != S2.c1 && return false
  change_base_ring(GF(2), A) * S1.w2 == S2.w2 || return false
  S2.p1 * A == S1.p1 || return false
  for i in 1:r, j in i:r, k in j:r
    _mu_eval(S2, A[:, i], A[:, j], A[:, k]) == _mu_entry(S1, i, j, k) || return false
  end
  return true
end

"""
    _transport_invariant_system(S, A)

Express `S` in new coordinates `y=A*x`, where `A` is integral unimodular.
This helper is useful for tests and for importing externally chosen bases.
"""
function _transport_invariant_system(S::_GKMInvariantSystem, A::ZZMatrix)
  r = S.H2_rank
  @req size(A) == (r, r) && abs(det(A)) == 1 "A must be an r by r unimodular matrix."
  Ainv_QQ = inv(change_base_ring(QQ, A))
  Ainv = zero_matrix(ZZ, r, r)
  for i in 1:r, j in 1:r
    @req denominator(Ainv_QQ[i, j]) == 1 "Inverse of a unimodular matrix must be integral."
    Ainv[i, j] = ZZ(numerator(Ainv_QQ[i, j]))
  end

  triples = _symmetric_triples(r)
  packed = zero_matrix(ZZ, 1, length(triples))
  for (q, (i, j, k)) in enumerate(triples)
    packed[1, q] = _mu_eval(S, Ainv[:, i], Ainv[:, j], Ainv[:, k])
  end
  Rz, z = polynomial_ring(ZZ, ["z$i" for i in 1:r])
  cubic = sum(
    _multiplicity(i, j, k) * packed[1, q] * z[i] * z[j] * z[k]
    for (q, (i, j, k)) in enumerate(triples); init=zero(Rz)
  )
  c1 = A * S.c1
  diagnostics = copy(S.diagnostics)
  transported = _GKMInvariantSystem(
    r,
    free_module(ZZ, r),
    S.basis_localizations * Ainv,
    packed,
    triples,
    Rz,
    cubic,
    c1,
    change_base_ring(GF(2), c1),
    S.p1 * Ainv,
    S.c2 * Ainv,
    S.b3,
    S.euler,
    S.characteristic_numbers,
    S.root,
    S.weight_rank,
    diagnostics,
    Dict{Symbol, Any}(),
  )
  @assert _validate_invariant_system(transported)
  @assert _verify_system_isomorphism(S, transported, A; preserve_almost_complex=true)
  return transported
end

function _content(M)
  c = ZZ(0)
  for x in M
    c = gcd(c, abs(ZZ(x)))
  end
  return c
end

function _vector_divisibility(v)
  return _content(v)
end

function _smith_diagonal(M::ZZMatrix)
  D = snf(M)
  result = ZZRingElem[]
  for i in 1:min(nrows(D), ncols(D))
    !iszero(D[i, i]) && push!(result, abs(D[i, i]))
  end
  return result
end

function _tensor_signature(S::_GKMInvariantSystem)
  haskey(S.cache, :tensor_signature) && return S.cache[:tensor_signature]
  F = _mu_flattening(S)
  result = (
    content=_content(S.mu_packed),
    rank=rank(change_base_ring(QQ, F)),
    radical_rank=S.H2_rank - rank(change_base_ring(QQ, F)),
    smith=_smith_diagonal(F),
  )
  S.cache[:tensor_signature] = result
  return result
end

function _rational_signature(B::ZZMatrix)
  # Symmetric Gaussian elimination by congruence. A nonzero diagonal pivot
  # contributes its sign; when every diagonal is zero, a nonzero off-diagonal
  # 2-by-2 pivot contributes one positive and one negative square.
  A = change_base_ring(QQ, B)
  positive = negative = zero_count = 0
  while nrows(A) > 0
    n = nrows(A)
    pivot = findfirst(i -> !iszero(A[i, i]), 1:n)
    if !isnothing(pivot)
      order = vcat([pivot], filter(!=(pivot), 1:n))
      A = A[order, order]
      a = A[1, 1]
      a > 0 ? (positive += 1) : (negative += 1)
      if n == 1
        A = zero_matrix(QQ, 0, 0)
      else
        schur = zero_matrix(QQ, n - 1, n - 1)
        for i in 1:(n - 1), j in 1:(n - 1)
          schur[i, j] = A[i + 1, j + 1] - A[i + 1, 1] * A[1, j + 1] / a
        end
        A = schur
      end
      continue
    end

    off_diagonal = nothing
    for i in 1:n, j in (i + 1):n
      if !iszero(A[i, j])
        off_diagonal = (i, j)
        break
      end
    end
    if isnothing(off_diagonal)
      zero_count += n
      break
    end

    i, j = off_diagonal
    order = vcat([i, j], filter(x -> x != i && x != j, 1:n))
    A = A[order, order]
    positive += 1
    negative += 1
    if n == 2
      A = zero_matrix(QQ, 0, 0)
    else
      b = A[1, 2]
      schur = zero_matrix(QQ, n - 2, n - 2)
      for i in 1:(n - 2), j in 1:(n - 2)
        correction = (A[1, i + 2] * A[2, j + 2] + A[2, i + 2] * A[1, j + 2]) / b
        schur[i, j] = A[i + 2, j + 2] - correction
      end
      A = schur
    end
  end
  return (positive, zero_count, negative)
end

function _contraction_signature(S::_GKMInvariantSystem, x)
  cache = get!(S.cache, :contraction_signature, Dict{Any, Any}())
  key = Tuple(ZZ(_vector_entry(x, i)) for i in 1:S.H2_rank)
  haskey(cache, key) && return cache[key]
  B = _contraction_matrix(S, x)
  rk = rank(change_base_ring(QQ, B))
  determinant = nrows(B) == 0 ? ZZ(1) : det(B)
  inertia = _rational_signature(B)
  result = (
    rank=rk,
    abs_determinant=abs(determinant),
    smith=_smith_diagonal(B),
    signature=inertia,
    even=all(i -> iseven(B[i, i]), 1:nrows(B)),
  )
  cache[key] = result
  return result
end

function _restricted_tensor_signature(S::_GKMInvariantSystem, K::ZZMatrix)
  s = ncols(K)
  triples = _symmetric_triples(s)
  packed = zero_matrix(ZZ, 1, length(triples))
  for (q, (i, j, k)) in enumerate(triples)
    packed[1, q] = _mu_eval(S, K[:, i], K[:, j], K[:, k])
  end
  pairs = [(j, k) for j in 1:s for k in j:s]
  flat = zero_matrix(ZZ, s, length(pairs))
  for i in 1:s, (q, (j, k)) in enumerate(pairs)
    t = _triple_index_from_list(triples, i, j, k)
    flat[i, q] = packed[1, t]
  end
  return (
    rank=s,
    content=_content(packed),
    flattening_rank=rank(change_base_ring(QQ, flat)),
    radical_rank=s - rank(change_base_ring(QQ, flat)),
    smith=_smith_diagonal(flat),
  )
end

function _canonical_sublattice_signatures(S::_GKMInvariantSystem; preserve_almost_complex::Bool)
  cache = get!(S.cache, :canonical_sublattice, Dict{Bool, Any}())
  haskey(cache, preserve_almost_complex) && return cache[preserve_almost_complex]
  r = S.H2_rank
  covectors = Pair{Symbol, ZZMatrix}[:p1 => S.p1]
  if preserve_almost_complex
    l11 = zero_matrix(ZZ, 1, r)
    for i in 1:r
      ei = zero_matrix(ZZ, r, 1)
      ei[i, 1] = 1
      l11[1, i] = _mu_eval(S, S.c1, S.c1, ei)
    end
    push!(covectors, :c1_squared => l11)
    push!(covectors, :c2 => S.c2)
  end

  result = Dict{Any, Any}()
  for (name, covector) in covectors
    K = _integral_kernel_basis(covector)
    ncols(K) >= 1 && (result[name] = _restricted_tensor_signature(S, K))
  end
  for i in 1:length(covectors), j in (i + 1):length(covectors)
    j > length(covectors) && continue
    name = (covectors[i].first, covectors[j].first)
    equations = vcat(covectors[i].second, covectors[j].second)
    K = _integral_kernel_basis(equations)
    ncols(K) >= 1 && (result[name] = _restricted_tensor_signature(S, K))
  end
  cache[preserve_almost_complex] = result
  return result
end

function _cheap_obstruction(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem,
  ; preserve_almost_complex::Bool,
)
  S1.H2_rank == S2.H2_rank || return "different H2 ranks"
  S1.b3 == S2.b3 || return "different b3"
  S1.euler == S2.euler || return "different Euler characteristics"
  _content(S1.p1) == _content(S2.p1) || return "different p1 contents"
  iszero(S1.w2) == iszero(S2.w2) || return "one w2 vanishes and the other does not"
  _tensor_signature(S1) == _tensor_signature(S2) || return "different tensor content, rank, radical, or flattening Smith data"
  _canonical_sublattice_signatures(S1; preserve_almost_complex) ==
    _canonical_sublattice_signatures(S2; preserve_almost_complex) ||
      return "different tensor data on a canonical saturated kernel"

  if preserve_almost_complex
    _vector_divisibility(S1.c1) == _vector_divisibility(S2.c1) || return "different c1 divisibilities"
    (S1.p1 * S1.c1)[1, 1] == (S2.p1 * S2.c1)[1, 1] || return "different p1(c1)"
    S1.characteristic_numbers.c1_cubed == S2.characteristic_numbers.c1_cubed || return "different c1^3"
    S1.characteristic_numbers.c1_c2 == S2.characteristic_numbers.c1_c2 || return "different c1*c2"
    _contraction_signature(S1, S1.c1) == _contraction_signature(S2, S2.c1) || return "different c1-contraction lattices"
    d1 = _vector_divisibility(S1.c1)
    d2 = _vector_divisibility(S2.c1)
    if !iszero(d1) && !iszero(d2)
      primitive1 = divexact(S1.c1, d1)
      primitive2 = divexact(S2.c1, d2)
      _contraction_signature(S1, primitive1) == _contraction_signature(S2, primitive2) ||
        return "different primitive-c1 contraction lattices"
    end
  end
  return nothing
end

function _partial_tensor_matches(S1, S2, columns)
  jmax = length(columns)
  for i in 1:jmax, j in i:jmax, k in j:jmax
    k == jmax || continue
    _mu_eval(S2, columns[i], columns[j], columns[k]) == _mu_entry(S1, i, j, k) || return false
  end
  return true
end

function _column_matrix(columns, ring=ZZ)
  isempty(columns) && return zero_matrix(ring, 0, 0)
  r = length(columns[1])
  A = zero_matrix(ring, r, length(columns))
  for j in 1:length(columns), i in 1:r
    A[i, j] = ring(_vector_entry(columns[j], i))
  end
  return A
end

function _bounded_integral_witness_search(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem,
  bounds;
  preserve_almost_complex::Bool,
  node_cap::Int=5_000_000,
  residue_classes::Vector{ZZMatrix}=ZZMatrix[],
  residue_modulus::ZZRingElem=ZZ(1),
)
  r = S1.H2_rank
  r == 0 && return zero_matrix(ZZ, 0, 0)
  ordered_bounds = sort!(unique(Int.(collect(bounds))))
  filter!(>(0), ordered_bounds)
  previous_bound = 0
  Q1 = preserve_almost_complex ? _contraction_matrix(S1, S1.c1) : nothing
  Q2 = preserve_almost_complex ? _contraction_matrix(S2, S2.c1) : nothing

  for bound in ordered_bounds
    nodes = Ref(0)
    columns = Vector{Vector{ZZRingElem}}()

    function search_column(j, compatible_classes, uses_new_shell)
      nodes[] >= node_cap && return nothing
      for tuple in Iterators.product(ntuple(_ -> (-bound):bound, r)...)
        nodes[] += 1
        nodes[] > node_cap && return nothing
        all(iszero, tuple) && continue
        foldl(gcd, (abs(value) for value in tuple); init=0) == 1 || continue
        y = ZZ.(collect(tuple))
        sum(S2.p1[1, i] * y[i] for i in 1:r; init=ZZ(0)) == S1.p1[1, j] || continue

        next_compatible = compatible_classes
        if !isempty(residue_classes)
          next_compatible = Int[]
          for class_index in compatible_classes
            residue = residue_classes[class_index]
            all(i -> mod(y[i] - residue[i, j], residue_modulus) == 0, 1:r) &&
              push!(next_compatible, class_index)
          end
          isempty(next_compatible) && continue
        end

        push!(columns, y)
        A_partial = _column_matrix(columns)
        rank(change_base_ring(QQ, A_partial)) == j || (pop!(columns); continue)
        _partial_tensor_matches(S1, S2, columns) || (pop!(columns); continue)
        if preserve_almost_complex
          contraction_matches = all(
            a -> sum(
              columns[a[1]][u] * Q2[u, v] * columns[a[2]][v]
              for u in 1:r for v in 1:r; init=ZZ(0),
            ) == Q1[a[1], a[2]],
            ((i, j) for i in 1:j),
          )
          contraction_matches || (pop!(columns); continue)
        end

        next_uses_new_shell = uses_new_shell || maximum(abs, tuple) > previous_bound

        if j == r
          A = A_partial
          if next_uses_new_shell && abs(det(A)) == 1 && _verify_system_isomorphism(
            S1, S2, A;
            preserve_almost_complex,
          )
            return A
          end
        else
          found = search_column(j + 1, next_compatible, next_uses_new_shell)
          !isnothing(found) && return found
        end
        pop!(columns)
      end
      return nothing
    end

    initial_classes = isempty(residue_classes) ? Int[] : collect(eachindex(residue_classes))
    witness = search_column(1, initial_classes, false)
    !isnothing(witness) && return witness
    previous_bound = bound
  end
  return nothing
end

function _field_int(x)
  return ZZ(lift(ZZ, x))
end

function _finite_field_histogram(
  S::_GKMInvariantSystem,
  p::Int;
  preserve_almost_complex::Bool,
)
  F = GF(p)
  r = S.H2_rank
  histogram = Dict{Any, Int}()
  c1 = [F(S.c1[i, 1]) for i in 1:r]
  w2 = p == 2 ? [F(lift(ZZ, S.w2[i, 1])) for i in 1:r] : elem_type(F)[]
  for tuple in Iterators.product(ntuple(_ -> 0:(p - 1), r)...)
    x = [F(a) for a in tuple]
    p1x = sum(F(S.p1[1, i]) * x[i] for i in 1:r; init=zero(F))
    muxxx = _mu_eval(S, x, x, x)
    B = _contraction_matrix(S, x; base_ring=F)
    rank_B = rank(B)
    det_B = r == 0 ? one(F) : det(B)
    key_values = Any[
      _field_int(p1x),
      _field_int(muxxx),
      rank_B,
      _field_int(det_B),
    ]
    if preserve_almost_complex
      push!(key_values, _field_int(_mu_eval(S, c1, x, x)))
      push!(key_values, _field_int(_mu_eval(S, c1, c1, x)))
    end
    if p == 2
      push!(key_values, _field_int(_mu_eval(S, w2, x, x)))
      push!(key_values, _field_int(_mu_eval(S, w2, w2, x)))
    end
    key = Tuple(key_values)
    histogram[key] = get(histogram, key, 0) + 1
  end
  return histogram
end

function _finite_field_signature_obstruction(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem,
  p::Int,
  ;
  preserve_almost_complex::Bool,
  point_cap::Int,
)
  F = GF(p)
  rank(_mu_flattening(S1; base_ring=F)) == rank(_mu_flattening(S2; base_ring=F)) ||
    return "tensor flattening ranks differ modulo $p"
  if preserve_almost_complex
    rank(_contraction_matrix(S1, S1.c1; base_ring=F)) ==
      rank(_contraction_matrix(S2, S2.c1; base_ring=F)) ||
        return "c1-contraction ranks differ modulo $p"
  end
  if ZZ(p)^S1.H2_rank <= point_cap
    _finite_field_histogram(S1, p; preserve_almost_complex) ==
      _finite_field_histogram(S2, p; preserve_almost_complex) ||
        return "decorated point histograms differ modulo $p"
  end
  return nothing
end

function _finite_field_isomorphism_search(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem,
  p::Int,
  ;
  preserve_almost_complex::Bool,
  node_cap::Int,
)
  F = GF(p)
  r = S1.H2_rank
  r == 0 && return (:found, [zero_matrix(F, 0, 0)], 1)
  columns = Vector{Vector{elem_type(F)}}()
  solutions = typeof(zero_matrix(F, r, r))[]
  nodes = Ref(0)
  capped = Ref(false)

  function search_column(j)
    for tuple in Iterators.product(ntuple(_ -> 0:(p - 1), r)...)
      nodes[] += 1
      if nodes[] > node_cap
        capped[] = true
        return nothing
      end
      any(!iszero, tuple) || continue
      y = [F(a) for a in tuple]
      sum(F(S2.p1[1, i]) * y[i] for i in 1:r; init=zero(F)) == F(S1.p1[1, j]) || continue
      push!(columns, y)
      A_partial = _column_matrix(columns, F)
      rank(A_partial) == j || (pop!(columns); continue)
      _partial_tensor_matches(S1, S2, columns) || (pop!(columns); continue)
      if j == r
        A = A_partial
        c1_ok = !preserve_almost_complex || A * change_base_ring(F, S1.c1) == change_base_ring(F, S2.c1)
        w2_ok = p != 2 || A * change_base_ring(F, S1.w2) == change_base_ring(F, S2.w2)
        if c1_ok && w2_ok && !iszero(det(A))
          push!(solutions, A)
          pop!(columns)
          return A
        end
      else
        found = search_column(j + 1)
        !isnothing(found) && return found
        capped[] && return nothing
      end
      isempty(columns) || pop!(columns)
      capped[] && return nothing
    end
    return nothing
  end

  witness = search_column(1)
  !isnothing(witness) && return (:found, solutions, nodes[])
  capped[] && return (:capped, solutions, nodes[])
  return (:none, solutions, nodes[])
end

function _lift_residue_matrix(A, p::Int)
  result = zero_matrix(ZZ, nrows(A), ncols(A))
  for i in 1:nrows(A), j in 1:ncols(A)
    result[i, j] = mod(ZZ(lift(ZZ, A[i, j])), p)
  end
  return result
end

function _crt_combine_residues(a::ZZRingElem, modulus::ZZRingElem, b::ZZRingElem, p::Int)
  modulus_mod_p = Int(mod(modulus, p))
  correction = mod(Int(mod(b - a, p)) * invmod(modulus_mod_p, p), p)
  return mod(a + modulus * correction, modulus * p)
end

function _crt_residue_classes(modular_data; class_cap::Int=100_000)
  usable_data = [data for data in modular_data if !isempty(data.solutions)]
  isempty(usable_data) && return (ZZMatrix[], ZZ(1), false)
  r = nrows(first(first(usable_data).solutions))
  classes = ZZMatrix[zero_matrix(ZZ, r, r)]
  modulus = ZZ(1)
  complete = true
  for data in usable_data
    next_classes = ZZMatrix[]
    for residue in classes, solution in data.solutions
      lifted = _lift_residue_matrix(solution, data.prime)
      combined = zero_matrix(ZZ, r, r)
      for i in 1:r, j in 1:r
        combined[i, j] = _crt_combine_residues(
          residue[i, j], modulus, lifted[i, j], data.prime,
        )
      end
      push!(next_classes, combined)
      if length(next_classes) >= class_cap
        complete = false
        break
      end
    end
    classes = next_classes
    modulus *= data.prime
    complete &= data.complete
    isempty(classes) && break
  end
  return (classes, modulus, complete)
end

function _qq_matrix_to_integral(A::QQMatrix)
  B = zero_matrix(ZZ, nrows(A), ncols(A))
  for i in 1:nrows(A), j in 1:ncols(A)
    denominator(A[i, j]) == 1 || return nothing
    B[i, j] = ZZ(numerator(A[i, j]))
  end
  return B
end

function _inverse_unimodular(A::ZZMatrix)
  inverse = _qq_matrix_to_integral(inv(change_base_ring(QQ, A)))
  @assert !isnothing(inverse)
  return inverse
end

function _reduced_search_system(S::_GKMInvariantSystem; preserve_almost_complex::Bool)
  r = S.H2_rank
  r == 0 && return (S, identity_matrix(ZZ, 0))
  source = if preserve_almost_complex
    Q = _contraction_matrix(S, S.c1)
    Q * transpose(Q)
  else
    F = _mu_flattening(S)
    F * transpose(F)
  end
  majorant = source + identity_matrix(ZZ, r)
  reduced_gram, row_transform = lll_gram_with_transform(majorant)
  @assert row_transform * majorant * transpose(row_transform) == reduced_gram
  @assert abs(det(row_transform)) == 1
  coordinate_map = _inverse_unimodular(transpose(row_transform))
  reduced = _transport_invariant_system(S, coordinate_map)
  return (reduced, coordinate_map)
end

function _definite_contraction_search(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem;
  preserve_almost_complex::Bool,
  group_order_cap::Int=100_000,
)
  preserve_almost_complex || return (:skipped, nothing, "c1 is not distinguished")
  Q1 = _contraction_matrix(S1, S1.c1)
  Q2 = _contraction_matrix(S2, S2.c1)
  (iszero(det(Q1)) || iszero(det(Q2))) && return (:skipped, nothing, "c1 contraction is singular")
  sig1, sig2 = _rational_signature(Q1), _rational_signature(Q2)
  sig1 == sig2 || return (:none, nothing, "definite contractions have different signatures")
  r = S1.H2_rank
  (sig1[1] == r || sig1[3] == r) || return (:skipped, nothing, "c1 contraction is not definite")
  L1 = integer_lattice(; gram=change_base_ring(QQ, Q1))
  L2 = integer_lattice(; gram=change_base_ring(QQ, Q2))
  @assert is_definite(L1) && is_definite(L2)

  isometric, T0_QQ = is_isometric_with_isometry(L1, L2)
  isometric || return (:none, nothing, "definite contraction lattices are not isometric")
  T0 = _qq_matrix_to_integral(T0_QQ)
  isnothing(T0) && return (:skipped, nothing, "OSCAR returned a nonintegral lattice isometry")

  order = automorphism_group_order(L1)
  order > group_order_cap && return (:skipped, nothing, "definite automorphism group order $order exceeds cap $group_order_cap")
  generators_QQ = automorphism_group_generators(L1; ambient_representation=false)
  generators = ZZMatrix[]
  for g in generators_QQ
    gz = _qq_matrix_to_integral(g)
    isnothing(gz) && return (:skipped, nothing, "OSCAR returned a nonintegral lattice automorphism")
    push!(generators, gz)
  end

  identity = identity_matrix(ZZ, r)
  group = ZZMatrix[identity]
  matrix_key(A) = Tuple(A[i, j] for i in 1:nrows(A) for j in 1:ncols(A))
  seen = Set([matrix_key(identity)])
  cursor = 1
  while cursor <= length(group)
    P = group[cursor]
    cursor += 1
    for g in generators
      Pg = P * g
      key = matrix_key(Pg)
      key in seen && continue
      push!(seen, key)
      push!(group, Pg)
      length(group) > group_order_cap && return (:skipped, nothing, "automorphism enumeration exceeded its cap")
    end
  end
  length(group) == Int(order) || return (:skipped, nothing, "automorphism generators did not enumerate the advertised complete group")

  # OSCAR uses row coordinates: T*Q2*T^t=Q1. Our convention is
  # A^t*Q2*A=Q1, hence A=transpose(T0)*transpose(P).
  for P in group
    A = transpose(T0) * transpose(P)
    if _verify_system_isomorphism(
      S1, S2, A;
      preserve_almost_complex,
    )
      return (:found, A, length(group))
    end
  end
  return (:none, nothing, "all $(length(group)) definite-lattice isometries were checked")
end

@doc raw"""
    compare_systems(S1, S2; ...)

Compare two systems of invariants as obtained from [`system_of_invariants_6d`](@ref).
By default this compares the oriented smooth Wall--Jupp--Žubr data
``(H^2,b_3,\mu,p_1,w_2)``. Set `preserve_almost_complex=true` to additionally
require preservation of ``c_1`` and hence of the encoded homotopy class of the
almost-complex structure. Consequently, a `:not_equivalent` result in the
almost-complex mode need not obstruct an oriented diffeomorphism after forgetting
the almost-complex structure.

There are three possible results:

- `:equivalent`: the two systems are isomorphic. This result always carries an exactly verified unimodular witness.
- `:not_equivalent`: the systems are not isomorphic in the selected comparison mode. This result always cites a rigorous obstruction.
- `:unknown`: the bounded search for isomorphisms of the given systems was unsuccessful.

# Optional arguments:
- `preserve_almost_complex::Bool`: `false` by default. If `true`, additionally require an isomorphism to map the first Chern class of `S1` to that of `S2`.
- `primes`: `[2, 3, 5, 7]` by default. Distinct primes used for rigorous finite-field obstructions and for congruence classes guiding the integral search.
- `finite_field_point_cap::Int`: `200_000` by default. Compute a complete decorated point histogram modulo `p` only if `p^H2_rank` does not exceed this cap.
- `finite_field_isomorphism_cap::Int`: `2_000_000` by default. Maximum nodes in each finite-field isomorphism search. The search stops at its first witness; complete exhaustion without one is rigorous, while reaching the cap is inconclusive.
- `integral_search_bounds`: `[1, 2, 3, 4]` by default. Successive entry bounds for the integral witness search. Exhausting them is inconclusive.
- `use_definite_contraction::Bool`: `true` by default. In almost-complex mode, completely enumerate isometries when the canonical `c1` contraction is definite.

# Example 1

Let us verify that $\mathbb{CP}^3$ and the full flag variety for $\mathbb{C}^3$ are not diffeomorphic.

```jldoctest
julia> G = projective_space(GKM_graph, 3);

julia> F = flag_variety(GKM_graph, [1,1,1]);

julia> SG = system_of_invariants_6d(G);

julia> SF = system_of_invariants_6d(F);

julia> compare_systems(SG, SF)
System comparison: not_equivalent (different H2 ranks)
```

The test fails because $G$ and $F$ have different second Betti numbers, which
are preserved by any diffeomorphism.

# Example 2
Let ``\Sigma_n`` be the Hirzebruch surface with parameter ``n`` and let
``X_n := \Sigma_n\times\mathbb{CP}^1``.
We will verify that $X_3$ and $X_4$ are not diffeomorphic while $X_3$ and $X_5$ are diffeomorphic.

```jldoctest compare_systems_hirzebruch
julia> P1 = projective_space(GKM_graph, 1);

julia> Hirzebruch = [gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, n)) for n in 1:5];

julia> Hirzebruch_times_P1 = [Hirzebruch[n] * P1 for n in 1:5];

julia> S3 = system_of_invariants_6d(Hirzebruch_times_P1[3]);

julia> S4 = system_of_invariants_6d(Hirzebruch_times_P1[4]);

julia> S5 = system_of_invariants_6d(Hirzebruch_times_P1[5]);

julia> compare_systems(S3, S4)
System comparison: not_equivalent (one w2 vanishes and the other does not)

julia> compare_systems(S3, S5)
System comparison: equivalent

julia> ans.witness
[-1    1   0]
[ 0   -1   0]
[ 0    0   1]
```

The second result carries a verified integral matrix in its `witness` field.

"""
function compare_systems(
  S1::_GKMInvariantSystem,
  S2::_GKMInvariantSystem;
  preserve_almost_complex::Bool=false,
  primes=[2, 3, 5, 7],
  finite_field_point_cap::Int=200_000,
  finite_field_isomorphism_cap::Int=2_000_000,
  integral_search_bounds=[1, 2, 3, 4],
  use_definite_contraction::Bool=true,
)
  _require_valid_invariant_system(S1, "the first")
  _require_valid_invariant_system(S2, "the second")
  all(p -> p > 1 && is_prime(p), primes) || throw(ArgumentError("primes must contain only primes"))
  length(unique(primes)) == length(primes) || throw(ArgumentError("primes must be distinct for CRT combination"))
  finite_field_point_cap >= 0 || throw(ArgumentError("finite_field_point_cap must be nonnegative"))
  finite_field_isomorphism_cap >= 0 || throw(ArgumentError("finite_field_isomorphism_cap must be nonnegative"))

  mode = preserve_almost_complex ? :oriented_smooth_almost_complex : :oriented_smooth
  diagnostics = Any[(:comparison_mode, mode)]
  obstruction = _cheap_obstruction(S1, S2; preserve_almost_complex)
  if !isnothing(obstruction)
    push!(diagnostics, (:cheap_obstruction, obstruction))
    return _SystemComparisonResult(:not_equivalent, nothing, obstruction, diagnostics)
  end

  try
    obstruction = _rational_polynomial_obstruction(S1, S2)
    if !isnothing(obstruction)
      certificate = (:rational_factorization, obstruction)
      push!(diagnostics, certificate)
      return _SystemComparisonResult(:not_equivalent, nothing, certificate, diagnostics)
    end
    push!(diagnostics, (:rational_factorization, :passed))
  catch error
    push!(diagnostics, (:rational_factorization, :skipped, sprint(showerror, error)))
  end

  for p in primes
    obstruction = _finite_field_signature_obstruction(
      S1, S2, p;
      preserve_almost_complex,
      point_cap=finite_field_point_cap,
    )
    if !isnothing(obstruction)
      certificate = (:finite_field_obstruction, p, obstruction)
      push!(diagnostics, certificate)
      return _SystemComparisonResult(:not_equivalent, nothing, certificate, diagnostics)
    end
  end

  reduced1, change1 = _reduced_search_system(S1; preserve_almost_complex)
  reduced2, change2 = _reduced_search_system(S2; preserve_almost_complex)
  push!(diagnostics, (:search_basis_reduction, change1, change2))

  modular_data = Any[]
  for p in primes
    status, solutions, nodes = _finite_field_isomorphism_search(
      reduced1, reduced2, p;
      preserve_almost_complex,
      node_cap=finite_field_isomorphism_cap,
    )
    push!(diagnostics, (:finite_field_search, p, status, nodes, length(solutions)))
    if status == :none
      certificate = (:finite_field_nonexistence, p, nodes)
      return _SystemComparisonResult(:not_equivalent, nothing, certificate, diagnostics)
    end
    # A first modular witness is enough to keep this prime viable. It supplies
    # only a heuristic CRT branch: we deliberately do not enumerate every
    # modular isomorphism after finding one.
    push!(modular_data, (prime=p, solutions=solutions, complete=false))
  end

  if use_definite_contraction
    status, witness, detail = _definite_contraction_search(
      reduced1, reduced2;
      preserve_almost_complex,
    )
    push!(diagnostics, (:definite_contraction, status, detail))
    if status == :found
      original_witness = _inverse_unimodular(change2) * witness * change1
      @assert _verify_system_isomorphism(
        S1, S2, original_witness;
        preserve_almost_complex,
      )
      return _SystemComparisonResult(:equivalent, original_witness, nothing, diagnostics)
    elseif status == :none
      certificate = (:definite_contraction, detail)
      return _SystemComparisonResult(:not_equivalent, nothing, certificate, diagnostics)
    end
  else
    push!(diagnostics, (:definite_contraction, :disabled))
  end

  # Possible future extension: a low-rank Gröbner unit-ideal test over QQ can
  # provide an additional rigorous non-existence certificate. It is deliberately
  # not part of the API until a predictable, tested implementation is available.

  residue_classes, residue_modulus, complete_residues = _crt_residue_classes(modular_data)
  push!(diagnostics, (
    :crt_residue_classes, length(residue_classes), residue_modulus,
    complete_residues ? :complete : :heuristic,
  ))
  witness = if isempty(residue_classes)
    _bounded_integral_witness_search(
      reduced1, reduced2, integral_search_bounds;
      preserve_almost_complex,
    )
  else
    _bounded_integral_witness_search(
      reduced1, reduced2, integral_search_bounds;
      preserve_almost_complex,
      residue_classes,
      residue_modulus,
    )
  end
  if isnothing(witness) && !complete_residues
    witness = _bounded_integral_witness_search(
      reduced1, reduced2, integral_search_bounds;
      preserve_almost_complex,
    )
  end
  if !isnothing(witness)
    original_witness = _inverse_unimodular(change2) * witness * change1
    @assert _verify_system_isomorphism(S1, S2, original_witness; preserve_almost_complex)
    push!(diagnostics, (:integral_search, :found, collect(integral_search_bounds)))
    return _SystemComparisonResult(:equivalent, original_witness, nothing, diagnostics)
  end
  push!(diagnostics, (
    :integral_search,
    :exhausted_bounds,
    collect(integral_search_bounds),
    :inconclusive,
    "the contraction was not a completely enumerable definite case and bounded integral search is not a non-existence proof",
  ))
  return _SystemComparisonResult(:unknown, nothing, nothing, diagnostics)
end
