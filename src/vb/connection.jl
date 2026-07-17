connection(E::AbstractGKMVectorBundle) = E.connection

"""
    build_vector_bundle_connection(G, GMtoM, weights) -> Connection

Build the natural connection on a GKM vector bundle. For a tangent bundle,
the transport is inherited from the chosen connection of its base.
"""
function build_vector_bundle_connection(
  G::AbstractGKMGraph{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}},
) where {R}
  bundle_transport = Dict{Edge,Vector{Int}}()
  bundle_coefficients = Dict{Edge,Vector{R}}()
  tangent_weights = _weights_are_tangent_weights(G, GMtoM, weights)

  for base_edge in edges(G), e in (base_edge, reverse(base_edge))
    image = tangent_weights ? copy(transport(connection(G))[e]) :
      _bundle_transport_along_edge(G, GMtoM, weights, e, R)
    bundle_transport[e] = image
    bundle_coefficients[e] = _bundle_coefficients(G, GMtoM, weights, e, image, R)
  end

  return Connection{R}(
    bundle_transport,
    bundle_coefficients,
    tangent_weights ? "Induced from base" : "Natural",
  )
end

function _weights_are_tangent_weights(G, GMtoM, weights)
  size(weights, 2) == valency(G) || return false
  return all(
    weights[v, i] == GMtoM(weight(G, v, i))
    for v in vertices(G) for i in 1:valency(G)
  )
end

function _bundle_transport_along_edge(G, GMtoM, weights, e::Edge, ::Type{R}) where {R}
  n = size(weights, 2)
  edge_weight = GMtoM(weight(G, e))
  candidates = Vector{Vector{Int}}(undef, n)
  for i in 1:n
    candidates[i] = [j for j in 1:n if !isnothing(_connection_coefficient(
      weights[src(e), i] - weights[dst(e), j], edge_weight, R,
    ))]
    isempty(candidates[i]) && throw(ArgumentError("no compatible fibre weight for summand $i along $e"))
  end
  image, used = zeros(Int, n), falses(n)
  _find_bundle_matching!(image, used, candidates, 1) ||
    throw(ArgumentError("the compatible fibre weights along $e do not form a bijection"))
  return image
end

function _find_bundle_matching!(image, used, candidates, i)
  i > length(candidates) && return true
  for j in candidates[i]
    used[j] && continue
    image[i], used[j] = j, true
    _find_bundle_matching!(image, used, candidates, i + 1) && return true
    used[j] = false
  end
  return false
end

function _bundle_coefficients(G, GMtoM, weights, e::Edge, image, ::Type{R}) where {R}
  edge_weight = GMtoM(weight(G, e))
  result = Vector{R}(undef, length(image))
  for i in eachindex(image)
    coefficient = _connection_coefficient(
      weights[src(e), i] - weights[dst(e), image[i]], edge_weight, R,
    )
    isnothing(coefficient) && throw(ArgumentError("invalid fibre transport for summand $i along $e"))
    result[i] = coefficient
  end
  return result
end
