
function _euler_class(G, v::Int, t::Vector{T}) where {T}
  res = zero(t[1])
  for i in 1:valency(G)
    res += _flag_weight_class(G, v, i, t)
  end
  return res
    
end

function _euler_class(G, v::Int)
  return _euler_class(G, v, gens_coeffRing(G))
end

function first_chern_class(G)
  val = valency(G)
  t = gens_coeffRing(G)
  e = gens_cohomRing(G)
  res = zero(e[1])

  for v in 1:n_vertices(G.core.g)
    localFactor = zero(t[1])
    # Iterate over all flags at this vertex (including standalone flags)
    for i in 1:val
      localFactor += _flag_weight_class(G.core, v, i, t)
    end
    res += localFactor * e[v]
  end
  return res
end

function point_class(G, v::Int)
  R = gens_cohomRing(G)
  return R[v]
end

function point_class(G, Vertexlabel::String)
  index = find_vertex_index(Vertexlabel, G)
  return point_class(G, index)
end

function _localize_at_vertex(G, c, e)
  success, quo = divides(c*e, e)
  @req success "class $c is invalid"
  return quo
end

function localize_at_vertex(G, c, v::Int)
  return _localize_at_vertex(G, c, gens_cohomRing(G)[v])
end

function localize_at_vertex(G, c, Vertexlabel::String)
  v = find_vertex_index(Vertexlabel, G)
  return _localize_at_vertex(G, c, gens_cohomRing(G)[v])
end


function _localized_vertex_coefficient(G, c, v::Int)
  localized = localize_at_vertex(G, c, v)
  representative = lift(localized)
  # Multiplication by the vertex idempotent removes every e-variable, so the
  # representative must be a constant polynomial over Frac(H_T^*(pt)).
  is_constant(representative) || throw(ArgumentError(
    "localization at vertex $v did not reduce to a coefficient",
  ))
  return constant_coefficient(representative)
end

_is_polynomial_fraction(f) = is_unit(denominator(f))

"""
    is_gkm_class(G, c) -> Bool

Return whether `c` belongs to the (unlocalized) GKM cohomology ring.  Every
vertex localization must lie in `H_T^*(pt)`, and for every edge `e = (v,w)`
the difference `c_v - c_w` must be divisible by `weight_class(G, e)` in that
polynomial ring.
"""
function is_gkm_class(G::AbstractGKMGraph, c)
  H = get_cohomology(G).localized_cohomology
  parent(c) == H || throw(ArgumentError(
    "the class does not belong to the localized cohomology ring of G",
  ))

  localizations = [
    _localized_vertex_coefficient(G, c, v)
    for v in vertices(G)
  ]
  all(_is_polynomial_fraction, localizations) || return false

  for e in edges(G)
    difference = localizations[src(e)] - localizations[dst(e)]
    quotient = difference / weight_class(G, e)
    _is_polynomial_fraction(quotient) || return false
  end
  return true
end

function integrate(G, c, e::Edge)
  v = src(e)
  w = dst(e)
  e_gens = gens_cohomRing(G)
  c_v = _localize_at_vertex(G, c, e_gens[v])
  c_w = _localize_at_vertex(G, c, e_gens[w])
  return (c_v - c_w) // weight_class(G, e)
end

function integrate(G, c, s1::String, s2::String)
  e = Edge(find_vertex_index(s1, G), find_vertex_index(s2, G))
  return integrate(G, c, e)
end

function integrate(G, c)
  e = gens_cohomRing(G)
  ans = zero(e[1])
  one_R = one(e[1])
  for v in num_vertices(G)
    euler = _euler_class(G, v)*one_R
    add!(ans, _localize_at_vertex(G, c, e[v]) // euler)
  end
  return ans
end
