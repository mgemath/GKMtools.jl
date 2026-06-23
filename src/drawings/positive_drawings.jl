###############################################################################
#
#   Positive Drawings of GKM Graphs
#
###############################################################################

@doc raw"""
    positive_drawing_representative(G::AbstractGKM_graph; strong::Bool=false) -> Union{Nothing, NamedTuple}

Determine whether `G` admits a *positive* drawing and, if so, return a generic
representative of it together with its admissibility and convexity analysis.

A drawing is *positive* if it is admissible and, for every edge `e`, the drawn
vector `pos(dst(e)) - pos(src(e))` is a **positive** rational multiple of the axial
weight `w_e`; equivalently `l_e(c) > 0` for every edge, where `l_e` is the linear
form on the drawing space `D(G)` computed by [`drawing_space`](@ref).

In the central hyperplane arrangement `H = ⋃_e ker(l_e)` inside `D(G)`, the positive
drawings form the single open chamber with the all-`+1` sign vector.  Since the sign
of `l_e` is independent of the orientation of `e` (reversing `e` negates both `w_e`
and the displacement), positivity is well-defined on the canonical orientation from
`edges(G.g)`.  This chamber is unique up to the global `C`/`-C` identification, so
there is at most one positive drawing representative.

If the positive chamber is empty, no positive drawing exists and `nothing` is
returned.  Otherwise a generic interior point of the chamber is chosen (the average
of the extreme rays of the closed positive cone) and converted to a representative.

# Arguments
- `strong::Bool=false`: if `false`, the `convex` field reports weak convexity
  ([`is_weakly_convex_drawing`](@ref)); if `true`, it reports strong convexity
  ([`is_strongly_convex_drawing`](@ref)).  Matches [`convex_drawing_representatives`](@ref).

# Returns
- `nothing` if no positive drawing exists (or, as a sanity safeguard, if the computed
  representative is not admissible — this should never happen).
- otherwise a `NamedTuple` `(positions, admissible, convex)` where
  - `positions :: Vector{Vector{QQFieldElem}}` is the representative (length
    `n_vertices(G.g)`, vertex 1 at the origin),
  - `admissible :: Bool` is the admissibility sanity check (always `true` here),
  - `convex :: Bool` is whether the representative is (weakly/strongly) convex.

Positivity and convexity are independent properties: a positive drawing need not be
convex, and a convex drawing need not be positive.

# Compatibility

To pass the representative to [`project_drawings`](@ref) or [`latex_drawing`](@ref),
wrap it in a one-element vector, e.g. `project_drawings([result.positions])`.

# Example
```jldoctest
julia> G = flag_variety(GKM_graph, [2, 1]);

julia> result = positive_drawing_representative(G);

julia> result.admissible
true

julia> length(result.positions)
3
```
"""
function positive_drawing_representative(G::AbstractGKM_graph; strong::Bool=false)
  n = n_vertices(G.g)
  d = rank_torus(G)

  basis, l_e_dict = drawing_space(G)
  r = ncols(basis)

  # No drawing freedom: every edge collapses, so no positive drawing exists.
  if r == 0
    println("No positive drawing exists: the drawing space D(G) is trivial.")
    return nothing
  end

  edge_list = collect(edges(G.g))
  L = [l_e_dict[e] for e in edge_list]

  # If some edge's linear form is identically zero, that edge collapses in every
  # drawing and no admissible (hence no positive) drawing exists.
  for le in L
    if all(c -> c == 0, le)
      println("No positive drawing exists: an edge collapses in every drawing.")
      return nothing
    end
  end

  # Interior point of the open positive chamber { c : l_e(c) > 0 for all e }.
  c = _positive_chamber_interior(L, r)
  if isnothing(c)
    println("No positive drawing exists: the positive chamber is empty.")
    return nothing
  end

  positions = _coords_to_positions(basis, c, n, d, r)

  # Admissibility sanity check: every edge displacement must be nonzero.
  # This always holds for a strict interior point of the positive chamber.
  admissible = all(edge_list) do e
    v, w = src(e), dst(e)
    any(k -> positions[w][k] != positions[v][k], 1:d)
  end
  if !admissible
    @warn "Positive chamber representative is not admissible; this should be impossible."
    return nothing
  end

  convex = strong ? is_strongly_convex_drawing(G, positions) : is_weakly_convex_drawing(G, positions)

  return (positions = positions, admissible = admissible, convex = convex)
end

"""
    _positive_chamber_interior(L::Vector{Vector{QQFieldElem}}, r::Int) -> Union{Nothing, Vector{QQFieldElem}}

Return an interior point of the open cone `{ x ∈ QQ^r : L_i · x > 0 for all i }`,
or `nothing` if that open cone is empty.

The closed cone `{ x : L_i · x ≥ 0 ∀ i }` is built with OSCAR's
`cone_from_inequalities` (which models `{x : A x ≤ 0}`, so the forms are negated).
The open cone is nonempty iff this closed cone is full-dimensional (`dim == r`), and
in that case the average of the extreme rays is a strictly-interior point.
"""
function _positive_chamber_interior(L::Vector{Vector{QQFieldElem}}, r::Int)
  n_hyp = length(L)
  # No edges: positivity is vacuous, use a generic point.
  n_hyp == 0 && return [QQ(1) for _ in 1:r]

  C = cone_from_inequalities(-L)   # { x : L_i · x ≥ 0 for all i }
  dim(C) < r && return nothing
  R = rays(C)
  isempty(R) && return nothing
  ray_vecs = [Vector{QQFieldElem}(rv) for rv in R]
  p = sum(ray_vecs) // QQ(length(ray_vecs))
  # Verify strict positivity (should always hold for an interior point).
  all(i -> sum(L[i][k] * p[k] for k in 1:r) > 0, 1:n_hyp) || return nothing
  return p
end
