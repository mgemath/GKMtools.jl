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
    we = [QQ(_w(G, e)[k]) for k in 1:d]
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
    we = [QQ(_w(G, e)[k]) for k in 1:d]
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
  L = [l_e_dict[e] for e in edge_list]

  # Evaluate linear form at point
  function eval_form(Lf::Vector{QQFieldElem}, p::Vector{QQFieldElem})
    isempty(Lf) && return QQ(0)
    return sum(Lf[k] * p[k] for k in 1:r)
  end

  qq_sign(x::QQFieldElem) = x > 0 ? 1 : (x < 0 ? -1 : 0)

  sign_vec(p::Vector{QQFieldElem}) = [qq_sign(eval_form(L[i], p)) for i in 1:m]

  # Canonical form: normalize so first nonzero entry is positive (identifies C and -C)
  function canonical_sv(sv::Vector{Int})
    idx = findfirst(s -> s != 0, sv)
    isnothing(idx) || sv[idx] > 0 ? sv : -sv
  end

  # If any l_e is identically zero, no admissible drawings exist
  for i in 1:m
    if all(c -> c == 0, L[i])
      return Vector{Vector{Vector{QQFieldElem}}}()
    end
  end

  # If no hyperplanes, the whole D(G) is one admissible component
  if m == 0
    # No hyperplanes: the whole D(G) is one admissible component
    pos = [[QQ(0) for _ in 1:d] for _ in 1:n]
    return [pos]
  end

  # Find a generic starting point in QQ^r not on any hyperplane
  function find_generic_point()
    # Try ±standard basis vectors
    for k in 1:r, s in [1, -1]
      p = fill(QQ(0), r)
      p[k] = QQ(s)
      all(x -> x != 0, sign_vec(p)) && return p
    end
    # Try small integer combinations (all sign patterns of basis vectors)
    r_eff = min(r, 10)
    for bits in 0:(2^r_eff - 1)
      p = fill(QQ(0), r)
      for k in 1:r_eff
        p[k] = isodd(bits >> (k - 1)) ? QQ(-1) : QQ(1)
      end
      all(x -> x != 0, sign_vec(p)) && return p
    end
    # Random integer fallback
    for _ in 1:1000
      p = [QQ(rand(-20:20)) for _ in 1:r]
      all(c -> c == 0, p) && continue
      all(x -> x != 0, sign_vec(p)) && return p
    end
    return nothing
  end

  # Given a point p in chamber σ, find a point in the adjacent chamber obtained
  # by crossing the hyperplane L[cross_idx] = 0.
  # We move along direction d = -σ_i * L[cross_idx] so that L[cross_idx] decreases to 0.
  function cross_hyperplane(p::Vector{QQFieldElem}, cross_idx::Int)
    Le = L[cross_idx]
    Lep = eval_form(Le, p)
    σi = qq_sign(Lep)
    σi == 0 && return nothing

    # d = -σi * Le; eval_form(Le, d) = -σi * ||Le||^2
    norm2 = sum(Le[k]^2 for k in 1:r)
    norm2 == 0 && return nothing
    Led = -σi * norm2          # eval_form(Le, d)
    t_star = -Lep // Led       # > 0: the crossing time

    # Find the nearest future hyperplane crossing at t > t_star
    min_future = nothing
    for j in 1:m
      j == cross_idx && continue
      Lf = L[j]
      Lfd = -σi * sum(Lf[k] * Le[k] for k in 1:r)  # eval_form(Lf, d)
      Lfd == 0 && continue
      Lfp = eval_form(Lf, p)
      t_f = -Lfp // Lfd
      t_f > t_star || continue
      if isnothing(min_future) || t_f < min_future
        min_future = t_f
      end
    end

    ε = isnothing(min_future) ? QQ(1) : (min_future - t_star) // 2
    t = t_star + ε
    return [p[k] + t * (-σi * Le[k]) for k in 1:r]
  end

  # BFS over chambers
  p0 = find_generic_point()
  isnothing(p0) && return Vector{Vector{Vector{QQFieldElem}}}()

  σ0 = sign_vec(p0)
  any(s -> s == 0, σ0) && return Vector{Vector{Vector{QQFieldElem}}}()

  chamber_reps = Dict{Vector{Int}, Vector{QQFieldElem}}()
  chamber_reps[canonical_sv(σ0)] = p0
  queue = [p0]

  while !isempty(queue)
    p = popfirst!(queue)
    for i in 1:m
      eval_form(L[i], p) == 0 && continue
      q = cross_hyperplane(p, i)
      isnothing(q) && continue
      σq = sign_vec(q)
      any(s -> s == 0, σq) && continue
      σq_canon = canonical_sv(σq)
      haskey(chamber_reps, σq_canon) && continue
      chamber_reps[σq_canon] = q
      push!(queue, q)
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
