###############################################################################
#
#   Admissible Drawings of GKM Graphs
#
###############################################################################

@doc raw"""
    drawing_space(G::AbstractGKM_graph) -> (QQMatrix, Dict{Edge, Vector{QQFieldElem}})

Compute the space `D(G)` of drawings of the GKM graph `G`, and the linear forms
`l_e` for each edge `e`.

A *drawing* of `G` is a map from the vertex set to `QQ^d` (where `d = rank_torus(G)`)
such that for each edge `e = Edge(v, w)`, the vector `pos(w) - pos(v)` is a rational
multiple of `G.w[e]` in `QQ^d`.  The space `D(G)` consists of drawings modulo
constant translation (vertex 1 is fixed at the origin).

# Returns

- `basis :: QQMatrix`: an `N × r` matrix whose columns form a QQ-basis for `D(G)`,
  where `N = (n_vertices - 1) * d` and `r = dim D(G)`.
  Entry `basis[(v-2)*d + k, i]` is the `k`-th coordinate of vertex `v ≥ 2`
  in the `i`-th basis drawing (vertex 1 is always at the origin).

- `l_e :: Dict{Edge, Vector{QQFieldElem}}`: for each edge `e`, the length-`r` vector
  representing the linear form `l_e : D(G) → QQ`.  If the drawing is given by
  coordinates `c ∈ QQ^r` (so vertex positions are `basis * c`), then
  `dot(l_e[e], c)` is the scalar `λ` with `pos(w) - pos(v) = λ · G.w[e]`.

# Example
```jldoctest
julia> G = flag_variety(GKM_graph, [2, 1]);

julia> basis, l_e = drawing_space(G);

julia> size(basis)
(4, 1)

julia> length(l_e)
3
```
"""
function drawing_space(G::AbstractGKM_graph)
  n = n_vertices(G.g)
  d = rank_torus(G)

  # Degenerate cases
  if n <= 1 || d == 0
    basis = zero_matrix(QQ, 0, 0)
    l_e_dict = Dict{Edge, Vector{QQFieldElem}}(e => QQFieldElem[] for e in edges(G.g))
    return basis, l_e_dict
  end

  N = (n - 1) * d   # number of free variables: positions of vertices 2..n, each in QQ^d

  # Variable index for component k of vertex v's position (v >= 2)
  var_idx(v::Int, k::Int) = (v - 2) * d + k

  edge_list = collect(edges(G.g))

  # Number of constraint rows: C(d,2) per edge
  nrows = length(edge_list) * d * (d - 1) ÷ 2

  # Build linear constraint matrix A such that D(G) = ker(A).
  # For each edge e = (v, w) with weight w_e in QQ^d:
  #   delta[k] = pos(w)[k] - pos(v)[k]  (0 for vertex 1)
  # Constraint for j < k:
  #   w_e[j] * delta[k] - w_e[k] * delta[j] = 0
  A = zero_matrix(QQ, nrows, N)
  row = 0
  for e in edge_list
    v, w = src(e), dst(e)
    we = [QQ(G.w[e][k]) for k in 1:d]
    for j in 1:d, k in j+1:d
      row += 1
      if w >= 2
        A[row, var_idx(w, j)] += we[k]
        A[row, var_idx(w, k)] -= we[j]
      end
      if v >= 2
        A[row, var_idx(v, j)] -= we[k]
        A[row, var_idx(v, k)] += we[j]
      end
    end
  end

  # Kernel of A = D(G)
  basis = if nrows == 0
    identity_matrix(QQ, N)
  else
    _, K = nullspace(A)
    K
  end

  r = ncols(basis)

  # Compute the linear form l_e for each edge.
  # l_e(drawing) = (pos(w)[k*] - pos(v)[k*]) / w_e[k*]
  # where k* is the first nonzero index of w_e.
  l_e_dict = Dict{Edge, Vector{QQFieldElem}}()
  for e in edge_list
    v, w = src(e), dst(e)
    we = [QQ(G.w[e][k]) for k in 1:d]
    kstar = findfirst(k -> we[k] != 0, 1:d)
    @assert !isnothing(kstar) "Edge $e has zero weight vector"
    le = Vector{QQFieldElem}(undef, r)
    for i in 1:r
      pos_w = (w >= 2) ? basis[var_idx(w, kstar), i] : QQ(0)
      pos_v = (v >= 2) ? basis[var_idx(v, kstar), i] : QQ(0)
      le[i] = (pos_w - pos_v) / we[kstar]
    end
    l_e_dict[e] = le
  end

  return basis, l_e_dict
end

@doc raw"""
    admissible_drawing_representatives(G::AbstractGKM_graph) -> Vector{Vector{Vector{QQFieldElem}}}

Enumerate one representative for each connected component of the space of admissible
drawings of `G`, up to overall sign.

An *admissible* drawing is one in which no two adjacent vertices coincide, i.e.,
`l_e(drawing) ≠ 0` for every edge `e`.  The space of admissible drawings is the
complement `D(G) \ H` of the hyperplane arrangement `H = ⋃_e ker(l_e)` in `D(G)`.
Since `H` is central (each hyperplane passes through the origin), components `C`
and `-C` are always distinct but are considered equivalent; only one is returned.

# Returns

A `Vector` of representatives.  Each representative is a
`Vector{Vector{QQFieldElem}}` of length `n_vertices(G.g)`, where entry `v` is a
vector of length `rank_torus(G)` giving the QQ-coordinates of vertex `v`.
Vertex 1 is always at the origin.

# Compatibility

To pass a representative to [`latex_drawing`](@ref) (which expects
`Vector{Tuple{Float64,Float64}}` for 2-dimensional drawings), convert via:

```julia
pos_float = [(Float64(p[1]), Float64(p[2])) for p in positions]
```

# Example
```jldoctest
julia> G = flag_variety(GKM_graph, [2, 1]);

julia> reps = admissible_drawing_representatives(G);

julia> length(reps)
1

julia> length(reps[1])
3
```
"""
function admissible_drawing_representatives(G::AbstractGKM_graph)
  n = n_vertices(G.g)
  d = rank_torus(G)

  basis, l_e_dict = drawing_space(G)
  r = ncols(basis)

  edge_list = collect(edges(G.g))
  m = length(edge_list)

  if r == 0
    return Vector{Vector{Vector{QQFieldElem}}}()
  end

  # l_e as vectors in QQ^r
  L_all = [l_e_dict[e] for e in edge_list]

  # If any l_e is identically zero, no admissible drawings exist
  for i in 1:m
    if all(c -> c == 0, L_all[i])
      return Vector{Vector{Vector{QQFieldElem}}}()
    end
  end

  # If no hyperplanes, the whole D(G) is one admissible component
  if m == 0
    pos = [[QQ(0) for _ in 1:d] for _ in 1:n]
    return [pos]
  end

  # Deduplicate proportional linear forms: two forms that are ANY nonzero scalar
  # multiple of each other define the same hyperplane and must share a sign-vector
  # entry.  Normalizing by the signed first nonzero entry (not its absolute value)
  # collapses both positive and negative multiples into one canonical form.
  function canonical_form(v::Vector{QQFieldElem})
    idx = findfirst(c -> c != 0, v)
    isnothing(idx) && return v
    return v .// v[idx]
  end

  seen = Dict{Vector{QQFieldElem}, Int}()
  L = Vector{Vector{QQFieldElem}}()
  for le in L_all
    cf = canonical_form(le)
    if !haskey(seen, cf)
      push!(L, le)
      seen[cf] = length(L)
    end
  end
  n_hyp = length(L)

  # Canonical form for sign vectors: normalize so first nonzero entry is positive (identifies C and -C)
  function canonical_sv(sv::Vector{Int})
    idx = findfirst(s -> s != 0, sv)
    isnothing(idx) || sv[idx] > 0 ? sv : -sv
  end

  # Given a sign vector σ ∈ {±1}^m, check whether the open cone
  #   { x ∈ QQ^r : σ_i * L_i(x) > 0  for all i }
  # is nonempty, and if so return an interior point.
  #
  # We use OSCAR's polyhedral cone { x : σ_i * L_i(x) ≥ 0 for all i }.
  # This closed cone is full-dimensional (dim == r) iff the open cone is nonempty
  # (for a central arrangement over QQ, the open cone is nonempty iff the closed
  # cone has full dimension).  An interior point is obtained as the sum of the
  # extreme rays, converted back to QQ.
  function find_chamber_interior(σ::Vector{Int})
    # Build inequality matrix: rows are σ_i * L_i  (cone = {x : A*x ≥ 0})
    # cone_from_inequalities(A) creates {x : Ax ≤ 0}; negate to get {x : σ_i L_i(x) ≥ 0}
    ineqs = [σ[i] .* L[i] for i in 1:n_hyp]
    C_cone = cone_from_inequalities(-ineqs)
    dim(C_cone) < r && return nothing
    # Interior point: average of extreme rays (each ray is a Vector{QQFieldElem})
    R = rays(C_cone)
    isempty(R) && return nothing
    ray_vecs = [Vector{QQFieldElem}(r_vec) for r_vec in R]
    p = sum(ray_vecs) // QQ(length(ray_vecs))
    # Verify (should always hold, but be safe)
    all(i -> sum(L[i][k] * p[k] for k in 1:r) * σ[i] > 0, 1:n_hyp) || return nothing
    return p
  end

  # BFS over sign vectors.  Starting sign vector: found by checking ±e_k and
  # small combinations until we land in a full-dimensional cone.
  function find_initial_sv()
    # Try ±standard basis vectors
    for k in 1:r, s in [1, -1]
      p = fill(QQ(0), r); p[k] = QQ(s)
      σ = [sum(L[i][j] * p[j] for j in 1:r) > 0 ? 1 : -1 for i in 1:n_hyp]
      !any(i -> sum(L[i][j] * p[j] for j in 1:r) == 0, 1:n_hyp) && return σ
    end
    # Try all ±1 combinations of basis vectors
    r_eff = min(r, 12)
    for bits in 0:(2^r_eff - 1)
      p = [isodd(bits >> (k - 1)) ? QQ(-1) : QQ(1) for k in 1:r_eff]
      append!(p, fill(QQ(1), r - r_eff))
      vals = [sum(L[i][k] * p[k] for k in 1:r) for i in 1:n_hyp]
      all(v -> v != 0, vals) && return [v > 0 ? 1 : -1 for v in vals]
    end
    return nothing
  end

  sv0 = find_initial_sv()
  isnothing(sv0) && return Vector{Vector{Vector{QQFieldElem}}}()

  chamber_reps = Dict{Vector{Int}, Vector{QQFieldElem}}()

  csν0 = canonical_sv(sv0)
  p0 = find_chamber_interior(sv0)
  isnothing(p0) && return Vector{Vector{Vector{QQFieldElem}}}()

  chamber_reps[csν0] = p0
  queue = [sv0]

  while !isempty(queue)
    σ = popfirst!(queue)
    for i in 1:n_hyp
      # Candidate neighbour: flip sign i
      σ_new = copy(σ); σ_new[i] = -σ[i]
      csν_new = canonical_sv(σ_new)
      haskey(chamber_reps, csν_new) && continue
      p_new = find_chamber_interior(σ_new)
      isnothing(p_new) && continue   # sign vector not realizable
      chamber_reps[csν_new] = p_new
      push!(queue, σ_new)
    end
  end

  # Convert basis coordinates to vertex positions in QQ^d
  function coords_to_positions(c::Vector{QQFieldElem})
    positions = Vector{Vector{QQFieldElem}}(undef, n)
    positions[1] = [QQ(0) for _ in 1:d]
    for v in 2:n
      positions[v] = [sum(basis[(v - 2) * d + k, i] * c[i] for i in 1:r) for k in 1:d]
    end
    return positions
  end

  result = Vector{Vector{Vector{QQFieldElem}}}()
  for (_, c) in chamber_reps
    push!(result, coords_to_positions(c))
  end
  return result
end

@doc raw"""
    project_drawings(reps::Vector{Vector{Vector{QQFieldElem}}}; projection::Union{Nothing, QQMatrix} = nothing) -> Vector{Vector{Tuple{Float64, Float64}}}

Project the output of [`admissible_drawing_representatives`](@ref) to 2D `Float64` coordinates
suitable for use with [`latex_drawing`](@ref).

Each drawing in `reps` gives vertex positions in `QQ^d`.  A `2 × d` linear
projection matrix is applied to each position, and the result is converted to
`Float64`.

# Arguments
- `reps`: output of `admissible_drawing_representatives(G)`.
- `projection`: an optional `2 × d` `QQMatrix`.  If `nothing` (the default), a
  random integer projection is chosen (identity for `d = 2`, identity-extended
  by zeros for `d < 2`, otherwise small random integers).

# Returns
A `Vector{Vector{Tuple{Float64,Float64}}}` of length `length(reps)`.
Entry `[c][v]` is the projected 2D position of vertex `v` in component `c`,
directly usable as the `positions` argument to `latex_drawing`.

# Example
```jldoctest
julia> G = flag_variety(GKM_graph, [2, 1]);

julia> reps = admissible_drawing_representatives(G);

julia> proj = project_drawings(reps);

julia> length(proj)
1

julia> length(proj[1])
3

julia> all(p -> p isa Tuple{Float64,Float64}, proj[1])
true
```
"""
function project_drawings(
  reps::Vector{Vector{Vector{QQFieldElem}}};
  projection::Union{Nothing, QQMatrix} = nothing
)::Vector{Vector{Tuple{Float64, Float64}}}

  isempty(reps) && return Vector{Vector{Tuple{Float64, Float64}}}()

  d = length(reps[1][1])

  # Build projection matrix P (2 × d)
  P = if !isnothing(projection)
    @req nrows(projection) == 2 && ncols(projection) == d "projection must be a 2×$d matrix"
    projection
  elseif d == 0
    zero_matrix(QQ, 2, 0)
  elseif d == 1
    P_d1 = zero_matrix(QQ, 2, 1)   # map x → (x, 0)
    P_d1[1, 1] = QQ(1)
    P_d1
  elseif d == 2
    identity_matrix(QQ, 2)
  else
    # Random small-integer projection
    P_rand = zero_matrix(QQ, 2, d)
    for i in 1:2, j in 1:d
      P_rand[i, j] = rand(-3:3)
    end
    P_rand
  end

  qq_to_float(x::QQFieldElem) = Float64(numerator(x)) / Float64(denominator(x))

  return map(reps) do positions
    map(positions) do pos_v
      x = sum(P[1, k] * pos_v[k] for k in 1:d; init = QQ(0))
      y = sum(P[2, k] * pos_v[k] for k in 1:d; init = QQ(0))
      (qq_to_float(x), qq_to_float(y))
    end
  end
end
