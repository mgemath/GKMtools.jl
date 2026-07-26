struct _ConvexBundleFunctor end
struct _ConcaveBundleFunctor end

@inline _decorated_graph(dt::GW_decorated_tree) = dt.tree
@inline _decorated_graph(dt::GW_decorated_graph) = dt.g
@inline _edge_degrees(dt::GW_decorated_tree, e::Edge) = (edgeMult(e, dt),)
@inline _edge_degrees(dt::GW_decorated_graph, e::Edge) = edgeMult(e, dt)
@inline _decorated_valency(dt::GW_decorated_tree, v::Int) = degree(dt.tree, v)
@inline _decorated_valency(dt::GW_decorated_graph, v::Int) = valency(v, dt)

@inline function _bundle_top_weight(E::AbstractGKMVectorBundle, v::Int, t)
  result = one(t[1])
  @inbounds for i in 1:rank(E)
    result *= _bundle_weight_class(E, v, i, t)
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a >= 0 "virtual_zero_section is only implemented for convex vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in 0:(Int(a) * d)
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_edge_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  G, t = dt.gkm, dt.context.t
  bundle_coefficients = coefficients(connection(E))
  bundle_rank = rank(E)
  for e in edges(_decorated_graph(dt))
    image_edge = imageOf(e, dt)
    source_vertex = src(image_edge)
    tangent_weight = _weight_class(G, image_edge, t)
    coefficients_along_edge = bundle_coefficients[image_edge]
    for d in _edge_degrees(dt, e)
      step = tangent_weight // d
      @inbounds for i in 1:bundle_rank
        a = coefficients_along_edge[i]
        @req a < 0 "derivated_functor is only implemented for concave vector bundles"
        fiber_weight = _bundle_weight_class(E, source_vertex, i, t)
        for k in (Int(a) * d + 1):-1
          result *= fiber_weight - k * step
        end
      end
    end
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConvexBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result = result // _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _bundle_vertex_contribution(result, dt, E::AbstractGKMVectorBundle, ::_ConcaveBundleFunctor)
  t = dt.context.t
  for v in vertices(_decorated_graph(dt))
    exponent = _decorated_valency(dt, v) - 1
    exponent <= 0 && continue
    result *= _bundle_top_weight(E, imageOf(v, dt), t)^exponent
  end
  return result
end

function _gw_bundle_functor(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle, mode)
  @req baseof(E) == dt.gkm "The vector bundle and decorated graph must have the same base"
  if dt isa GW_decorated_graph
    @req all(iszero, dt.genus) "Only implemented for genus zero"
  end
  result = _bundle_edge_contribution(one(dt.context.t[1]), dt, E, mode)
  return _bundle_vertex_contribution(result, dt, E, mode)
end

"""
    virtual_zero_section(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the top Chern class of the genus-zero pushforward of `E`.
The bundle must be convex along every GKM edge.
"""
function virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConvexBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConvexBundleFunctor()))
end

"""
    derivated_functor(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the top Chern class of the genus-zero first derived pushforward of `E`.
The bundle must be concave along every GKM edge.
"""
function derivated_functor(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_gw_bundle_functor(dt, $E, _ConcaveBundleFunctor()))
  return EquivariantClass(rule, dt -> _gw_bundle_functor(dt, E, _ConcaveBundleFunctor()))
end

"""
    reduced_virtual_zero_section(E::AbstractGKMVectorBundle) -> EquivariantClass

Return the virtual zero section with the top Chern class at the last marked point divided out.
"""
function reduced_virtual_zero_section(E::AbstractGKMVectorBundle)::EquivariantClass
  rule = :(_reduced_virtual_zero_section(dt, $E))
  return EquivariantClass(rule, dt -> _reduced_virtual_zero_section(dt, E))
end

function _reduced_virtual_zero_section(dt::Union{GW_decorated_tree,GW_decorated_graph}, E::AbstractGKMVectorBundle)
  @req !isempty(dt.marks) "Need at least one marked point to reduce the virtual zero section"
  result = _gw_bundle_functor(dt, E, _ConvexBundleFunctor())
  marked_vertex = imageOf(last(dt.marks), dt)
  return result // _bundle_top_weight(E, marked_vertex, dt.context.t)
end
