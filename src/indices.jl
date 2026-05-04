###############################################################################
#
# Index-increasing generic directions for GKM graphs
#
# Given a GKM graph G with weight lattice N = G.M, a vector xi in N is called
# _generic_ if for each flag f of G with weight alpha we have <xi, alpha> != 0,
# where <.,.> is the standard Euclidean inner product.
#
# For generic xi, each vertex has a _xi-index_ = number of flags at that vertex
# with weight pairing negatively against xi. Orient each edge e={v,w} as (v,w)
# so that <xi, w(v,w)> > 0. xi is _index increasing_ (resp. _weakly index
# increasing_) if index(w) > index(v) (resp. >=) for every edge so oriented.
#
###############################################################################

# Euclidean pairing of two weight-lattice elements, returned as QQFieldElem.
function _xi_pair(xi::AbstractAlgebra.Generic.FreeModuleElem{R},
                  alpha::AbstractAlgebra.Generic.FreeModuleElem{R})::QQFieldElem where R <: GKM_weight_type
  r = rank(parent(xi))
  s = QQ(0)
  for j in 1:r
    s += QQ(xi[j]) * QQ(alpha[j])
  end
  return s
end

@doc raw"""
    is_generic(G::AbstractGKM_graph, xi) -> Bool

Return `true` if `xi` is a generic element of the weight lattice `G.M`, i.e.
the Euclidean pairing of `xi` with every flag weight of `G` is nonzero.

# Example
```jldoctest is_generic_xi
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> g1, g2, g3 = gens(P2.M);

julia> is_generic(P2, g1)
false

julia> is_generic(P2, g1 + g2)
false

julia> is_generic(P2, g1 + 2*g2 + 4*g3)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> h1, h2, h3 = gens(P2_Q.M);

julia> is_generic(P2_Q, h1 + 2*h2 + 4*h3)
true
```
"""
function is_generic(G::AbstractGKM_graph{R},
                    xi::AbstractAlgebra.Generic.FreeModuleElem{R})::Bool where R <: GKM_weight_type
  @req parent(xi) === G.M "xi must live in G.M"
  for v in 1:n_vertices(G.g)
    for alpha in G.weights_at_vertex[v]
      if iszero(_xi_pair(xi, alpha))
        return false
      end
    end
  end
  return true
end

@doc raw"""
    xi_index(G::AbstractGKM_graph, xi, v::Int64) -> Int64

Return the `xi`-index of vertex `v` of `G`, i.e. the number of flags at `v`
whose weight pairs negatively with `xi`. Requires `xi` to be generic at `v`.

# Example
```jldoctest xi_index_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> g1, g2, g3 = gens(P2.M);

julia> xi = g1 + 2*g2 + 4*g3;

julia> [xi_index(P2, xi, v) for v in 1:3]
3-element Vector{Int64}:
 2
 1
 0

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> h1, h2, h3 = gens(P2_Q.M);

julia> [xi_index(P2_Q, h1 + 2*h2 + 4*h3, v) for v in 1:3]
3-element Vector{Int64}:
 2
 1
 0
```
"""
function xi_index(G::AbstractGKM_graph{R},
                  xi::AbstractAlgebra.Generic.FreeModuleElem{R},
                  v::Int64)::Int64 where R <: GKM_weight_type
  @req parent(xi) === G.M "xi must live in G.M"
  @req 1 <= v <= n_vertices(G.g) "Vertex $v out of bounds"
  cnt = 0
  for alpha in G.weights_at_vertex[v]
    p = _xi_pair(xi, alpha)
    @req !iszero(p) "xi is not generic at vertex $v"
    if p < 0
      cnt += 1
    end
  end
  return cnt
end

@doc raw"""
    is_index_increasing(G::AbstractGKM_graph, xi) -> Bool

Return `true` if `xi` is generic and index increasing on `G`.
For each edge `e = {v, w}`, orient it as `(v, w)` so that the weight of
`Edge(v, w)` pairs positively with `xi`; then require the `xi`-index at `w`
to be strictly greater than the `xi`-index at `v`.

# Example
```jldoctest is_index_increasing_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> g1, g2, g3 = gens(P2.M);

julia> is_index_increasing(P2, g1 + 2*g2 + 4*g3)
true

julia> is_index_increasing(P2, g1)
false

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> h1, h2, h3 = gens(P2_Q.M);

julia> is_index_increasing(P2_Q, h1 + 2*h2 + 4*h3)
true
```
"""
function is_index_increasing(G::AbstractGKM_graph{R},
                             xi::AbstractAlgebra.Generic.FreeModuleElem{R})::Bool where R <: GKM_weight_type
  return _check_increasing(G, xi, true)
end

@doc raw"""
    is_weakly_index_increasing(G::AbstractGKM_graph, xi) -> Bool

Return `true` if `xi` is generic and weakly index increasing on `G`.
For each edge `e = {v, w}`, orient it as `(v, w)` so that the weight of
`Edge(v, w)` pairs positively with `xi`; then require the `xi`-index at `w`
to be greater than or equal to the `xi`-index at `v`.

# Example
```jldoctest is_weakly_index_increasing_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> g1, g2, g3 = gens(P2.M);

julia> is_weakly_index_increasing(P2, g1 + 2*g2 + 4*g3)
true

julia> is_weakly_index_increasing(P2, g1)
false

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> h1, h2, h3 = gens(P2_Q.M);

julia> is_weakly_index_increasing(P2_Q, h1 + 2*h2 + 4*h3)
true
```
"""
function is_weakly_index_increasing(G::AbstractGKM_graph{R},
                                    xi::AbstractAlgebra.Generic.FreeModuleElem{R})::Bool where R <: GKM_weight_type
  return _check_increasing(G, xi, false)
end

function _check_increasing(G::AbstractGKM_graph{R},
                           xi::AbstractAlgebra.Generic.FreeModuleElem{R},
                           strict::Bool)::Bool where R <: GKM_weight_type
  @req parent(xi) === G.M "xi must live in G.M"
  is_generic(G, xi) || return false
  nv = n_vertices(G.g)
  indices = [xi_index(G, xi, v) for v in 1:nv]
  for e in edges(G.g)
    p = _xi_pair(xi, _w(G, e))
    @assert !iszero(p)
    v_from, v_to = p > 0 ? (src(e), dst(e)) : (dst(e), src(e))
    if strict
      indices[v_to] > indices[v_from] || return false
    else
      indices[v_to] >= indices[v_from] || return false
    end
  end
  return true
end

# Canonicalize a nonzero rational vector to a primitive integer vector whose
# first nonzero entry is positive. This identifies weights that differ by a
# nonzero rational scalar (same hyperplane).
function _canonical_normal(vec::Vector{QQFieldElem})::Vector{QQFieldElem}
  d = ZZ(1)
  for x in vec
    d = lcm(d, denominator(x))
  end
  ints = [numerator(x * d) for x in vec]
  g = ZZ(0)
  for x in ints
    g = gcd(g, x)
  end
  if !iszero(g)
    ints = [divexact(x, g) for x in ints]
  end
  for x in ints
    if !iszero(x)
      if x < 0
        ints = [-y for y in ints]
      end
      break
    end
  end
  return QQFieldElem[QQ(x) for x in ints]
end

function _has_zero_flag_weight(G::AbstractGKM_graph)::Bool
  for v in 1:n_vertices(G.g)
    for alpha in G.weights_at_vertex[v]
      iszero(alpha) && return true
    end
  end
  return false
end

function _collect_hyperplane_normals(G::AbstractGKM_graph)::Vector{Vector{QQFieldElem}}
  r = rank_torus(G)
  normals = Vector{Vector{QQFieldElem}}()
  seen = Set{Vector{QQFieldElem}}()
  for v in 1:n_vertices(G.g)
    for alpha in G.weights_at_vertex[v]
      vec = QQFieldElem[QQ(alpha[j]) for j in 1:r]
      all(iszero, vec) && continue
      canon = _canonical_normal(vec)
      if !(canon in seen)
        push!(seen, canon)
        push!(normals, canon)
      end
    end
  end
  return normals
end

function _chamber_cone(normals::Vector{Vector{QQFieldElem}}, s::Vector{Int8})
  nH = length(normals)
  r = length(normals[1])
  A = zero_matrix(QQ, nH, r)
  for i in 1:nH
    for j in 1:r
      A[i, j] = s[i] * normals[i][j]
    end
  end
  return cone_from_inequalities(A)
end

function _cone_representative(C, r::Int64)::Vector{QQFieldElem}
  rs = collect(rays(C))
  if isempty(rs)
    rs = rays_modulo_lineality(C).rays_modulo_lineality
  end
  out = QQFieldElem[QQ(0) for _ in 1:r]
  for ray in rs
    for i in 1:r
      out[i] += QQ(ray[i])
    end
  end
  return out
end

function _initial_sign_vector(normals::Vector{Vector{QQFieldElem}}, r::Int64)::Vector{Int8}
  nH = length(normals)
  # The moment curve meets each hyperplane in at most r - 1 points, so
  # checking nH * (r - 1) + 1 parameters guarantees a generic direction.
  for t in 0:(nH * (r - 1))
    xi = [ZZ(t)^(j - 1) for j in 1:r]
    s = Vector{Int8}(undef, nH)
    ok = true
    for i in 1:nH
      p = QQ(0)
      for j in 1:r
        p += xi[j] * normals[i][j]
      end
      if iszero(p)
        ok = false
        break
      end
      s[i] = p > 0 ? Int8(1) : Int8(-1)
    end
    ok && return s
  end
  error("Could not find a generic starting direction on the moment curve")
end

# Convert a QQ vector to a free-module element of G.M, clearing denominators
# and dividing out common factors so the result is a primitive element of the
# weight lattice. Falls through unchanged for QQ weight lattices.
function _vec_to_M(G::AbstractGKM_graph{R},
                   v::Vector{QQFieldElem})::AbstractAlgebra.Generic.FreeModuleElem{R} where R <: GKM_weight_type
  r = rank_torus(G)
  M = G.M
  if R == ZZRingElem
    d = ZZ(1)
    for x in v
      d = lcm(d, denominator(x))
    end
    ints = [numerator(x * d) for x in v]
    g = ZZ(0)
    for x in ints
      g = gcd(g, x)
    end
    if !iszero(g)
      ints = [divexact(x, g) for x in ints]
    end
    return sum(j -> ints[j] * gens(M)[j], 1:r)
  else
    return sum(j -> v[j] * gens(M)[j], 1:r)
  end
end

# BFS over chamber sign vectors. If `stop_on` is supplied, stop as soon as it
# returns `true` for some chamber representative and return just that one.
function _enumerate_chambers(G::AbstractGKM_graph{R},
                             normals::Vector{Vector{QQFieldElem}};
                             stop_on::Union{Nothing, Function} = nothing
                            )::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  r = rank_torus(G)
  if isempty(normals)
    return [gens(G.M)[1]]
  end
  init = _initial_sign_vector(normals, r)
  out = Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}()
  seen = Set{Vector{Int8}}()
  push!(seen, init)
  queue = Vector{Vector{Int8}}()
  push!(queue, init)
  while !isempty(queue)
    s = popfirst!(queue)
    C = _chamber_cone(normals, s)
    dim(C) < r && continue
    rep_vec = _cone_representative(C, r)
    rep = _vec_to_M(G, rep_vec)
    if !isnothing(stop_on) && stop_on(rep)
      return [rep]
    end
    if isnothing(stop_on)
      push!(out, rep)
    end
    for i in 1:length(s)
      s_new = copy(s)
      s_new[i] = -s_new[i]
      if !(s_new in seen)
        push!(seen, s_new)
        push!(queue, s_new)
      end
    end
  end
  return out
end

@doc raw"""
    generic_xi_representatives(G::AbstractGKM_graph) -> Vector

Return one representative per connected component of the set of generic
directions `xi` in the weight lattice of `G`. Each representative lives in
`G.M` (as a primitive element of the lattice if `G.M` is a free `ZZ`-module).

The components of the generic locus are the open chambers of the hyperplane
arrangement cut out by the flag weights of `G`, viewed in `G.M` tensored
with `QQ`.

# Example
On `P^1` there is a single hyperplane so two chambers:
```jldoctest generic_xi_representatives_example
julia> P1 = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> generic_xi_representatives(P1)
2-element Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}:
 (-1, 0)
 (0, -1)

julia> P1_Q = convert_weights(P1); # Let's also test with QQ-weights

julia> generic_xi_representatives(P1_Q)
2-element Vector{AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}:
 (-1, 0)
 (0, -1)
```
On `P^2` there are 6 chambers corresponding to the 6 orderings of the three
coordinates of `xi`:
```jldoctest generic_xi_representatives_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> reps = generic_xi_representatives(P2);

julia> length(reps)
6

julia> all(xi -> is_generic(P2, xi), reps)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> reps_Q = generic_xi_representatives(P2_Q);

julia> length(reps_Q)
6

julia> all(xi -> is_generic(P2_Q, xi), reps_Q)
true
```
Non-compact example: the total space of `O(1) + O(-1)` on `P^1`:
```jldoctest generic_xi_representatives_example
julia> T = total_space(vector_bundle_O(1, [1, -1]))
GKM graph with 2 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
Standalone flags:
1.2 => (0, 0, 1, 0)
1.3 => (0, 0, 0, 1)
2.2 => (-1, 1, 1, 0)
2.3 => (1, -1, 0, 1)

julia> length(generic_xi_representatives(T))
18

julia> T_Q = convert_weights(T); # Let's also test with QQ-weights

julia> length(generic_xi_representatives(T_Q))
18
```
"""
function generic_xi_representatives(G::AbstractGKM_graph{R}
    )::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  @req rank_torus(G) >= 1 "Torus has rank zero"
  _has_zero_flag_weight(G) && return AbstractAlgebra.Generic.FreeModuleElem{R}[]
  normals = _collect_hyperplane_normals(G)
  return _enumerate_chambers(G, normals)
end

@doc raw"""
    index_increasing_xi_representatives(G::AbstractGKM_graph) -> Vector

Return one representative per connected component of the generic locus on
which `xi` is index increasing. Each representative lives in `G.M`.

# Example
Every chamber of `P^2` is index-increasing:
```jldoctest index_increasing_xi_representatives_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> reps = index_increasing_xi_representatives(P2);

julia> length(reps)
6

julia> all(xi -> is_index_increasing(P2, xi), reps)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> reps_Q = index_increasing_xi_representatives(P2_Q);

julia> length(reps_Q)
6

julia> all(xi -> is_index_increasing(P2_Q, xi), reps_Q)
true
```
For the non-compact total space of `O(1) + O(-1)` on `P^1`, only some chambers
are strictly index-increasing (compare the count with the 18 generic chambers):
```jldoctest index_increasing_xi_representatives_example
julia> T = total_space(vector_bundle_O(1, [1, -1]));

julia> length(index_increasing_xi_representatives(T))
14

julia> T_Q = convert_weights(T); # Let's also test with QQ-weights

julia> length(index_increasing_xi_representatives(T_Q))
14
```
"""
function index_increasing_xi_representatives(G::AbstractGKM_graph{R}
    )::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  return filter(xi -> is_index_increasing(G, xi), generic_xi_representatives(G))
end

@doc raw"""
    weakly_index_increasing_xi_representatives(G::AbstractGKM_graph) -> Vector

Return one representative per connected component of the generic locus on
which `xi` is weakly index increasing. Each representative lives in `G.M`.

# Example
```jldoctest weakly_index_increasing_xi_representatives_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> reps = weakly_index_increasing_xi_representatives(P2);

julia> length(reps)
6

julia> all(xi -> is_weakly_index_increasing(P2, xi), reps)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> reps_Q = weakly_index_increasing_xi_representatives(P2_Q);

julia> length(reps_Q)
6

julia> all(xi -> is_weakly_index_increasing(P2_Q, xi), reps_Q)
true
```
On the non-compact total space of `O(1) + O(-1)` on `P^1`, every generic
chamber turns out to be weakly index-increasing even though only 14 of them
are strictly index-increasing:
```jldoctest weakly_index_increasing_xi_representatives_example
julia> T = total_space(vector_bundle_O(1, [1, -1]));

julia> length(weakly_index_increasing_xi_representatives(T))
18

julia> T_Q = convert_weights(T); # Let's also test with QQ-weights

julia> length(weakly_index_increasing_xi_representatives(T_Q))
18
```
"""
function weakly_index_increasing_xi_representatives(G::AbstractGKM_graph{R}
    )::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  return filter(xi -> is_weakly_index_increasing(G, xi), generic_xi_representatives(G))
end

@doc raw"""
    admits_index_increasing_xi(G::AbstractGKM_graph) -> Tuple{Bool, FreeModuleElem}

Return `(true, xi)` if `G` admits some generic index-increasing direction
`xi` in `G.M`, where `xi` is the first such direction found during chamber
enumeration. Return `(false, zero(G.M))` otherwise.

# Example
```jldoctest admits_index_increasing_xi_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> ok, xi = admits_index_increasing_xi(P2);

julia> ok
true

julia> is_index_increasing(P2, xi)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> ok_Q, xi_Q = admits_index_increasing_xi(P2_Q);

julia> ok_Q && is_index_increasing(P2_Q, xi_Q)
true

julia> G = gkm_2d([1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1]) # Blowup of P1 x P1 in 2 points
GKM graph with 6 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (-1, -1)
4 -> 3 => (0, -1)
5 -> 4 => (1, 0)
6 -> 1 => (0, -1)
6 -> 5 => (1, 1)

julia> admits_index_increasing_xi(G)
(false, (0, 0))
```
The same is true on the non-compact total space of `O(1) + O(-1)` on `P^1`:
```jldoctest admits_index_increasing_xi_example
julia> T = total_space(vector_bundle_O(1, [1, -1]));

julia> ok, xi = admits_index_increasing_xi(T);

julia> ok
true

julia> is_index_increasing(T, xi)
true

julia> T_Q = convert_weights(T); # Let's also test with QQ-weights

julia> ok_Q, xi_Q = admits_index_increasing_xi(T_Q);

julia> ok_Q && is_index_increasing(T_Q, xi_Q)
true
```
"""
function admits_index_increasing_xi(G::AbstractGKM_graph{R}
    )::Tuple{Bool, AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  @req rank_torus(G) >= 1 "Torus has rank zero"
  _has_zero_flag_weight(G) && return (false, zero(G.M))
  normals = _collect_hyperplane_normals(G)
  res = _enumerate_chambers(G, normals;
                            stop_on = xi -> is_index_increasing(G, xi))
  isempty(res) && return (false, zero(G.M))
  return (true, res[1])
end

@doc raw"""
    admits_weakly_index_increasing_xi(G::AbstractGKM_graph) -> Tuple{Bool, FreeModuleElem}

Return `(true, xi)` if `G` admits some generic weakly-index-increasing
direction `xi` in `G.M`, where `xi` is the first such direction found during
chamber enumeration. Return `(false, zero(G.M))` otherwise.

# Example
```jldoctest admits_weakly_index_increasing_xi_example
julia> P2 = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> ok, xi = admits_weakly_index_increasing_xi(P2);

julia> ok
true

julia> is_weakly_index_increasing(P2, xi)
true

julia> P2_Q = convert_weights(P2); # Let's also test with QQ-weights

julia> ok_Q, xi_Q = admits_weakly_index_increasing_xi(P2_Q);

julia> ok_Q && is_weakly_index_increasing(P2_Q, xi_Q)
true

julia> G = gkm_2d([1 0; 1 1; 0 1; -1 0; -1 -1; 0 -1]) # Blowup of P1 x P1 in 2 points
GKM graph with 6 nodes, valency 2 and axial function:
2 -> 1 => (-1, 0)
3 -> 2 => (-1, -1)
4 -> 3 => (0, -1)
5 -> 4 => (1, 0)
6 -> 1 => (0, -1)
6 -> 5 => (1, 1)

julia> admits_weakly_index_increasing_xi(G)
(true, (-1, -1))
```
Non-compact total space of `O(1) + O(-1)` on `P^1`:
```jldoctest admits_weakly_index_increasing_xi_example
julia> T = total_space(vector_bundle_O(1, [1, -1]));

julia> ok, xi = admits_weakly_index_increasing_xi(T);

julia> ok
true

julia> is_weakly_index_increasing(T, xi)
true

julia> T_Q = convert_weights(T); # Let's also test with QQ-weights

julia> ok_Q, xi_Q = admits_weakly_index_increasing_xi(T_Q);

julia> ok_Q && is_weakly_index_increasing(T_Q, xi_Q)
true
```
"""
function admits_weakly_index_increasing_xi(G::AbstractGKM_graph{R}
    )::Tuple{Bool, AbstractAlgebra.Generic.FreeModuleElem{R}} where R <: GKM_weight_type
  @req rank_torus(G) >= 1 "Torus has rank zero"
  _has_zero_flag_weight(G) && return (false, zero(G.M))
  normals = _collect_hyperplane_normals(G)
  res = _enumerate_chambers(G, normals;
                            stop_on = xi -> is_weakly_index_increasing(G, xi))
  isempty(res) && return (false, zero(G.M))
  return (true, res[1])
end
