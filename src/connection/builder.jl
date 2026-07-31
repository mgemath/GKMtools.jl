@doc raw"""
    build_gkm_connection(data::GKMCombinatorialData) -> Connection

Build an algorithmic connection compatible with the axial function of a
smooth or orbifold GKM graph. Along `e = (v,w)` it satisfies
`weight(G,w,transport[e][i]) = weight(G,v,i) - a[e][i]*weight(G,e)`.

If several transports are possible, return the lexicographically first one.
Throw `ArgumentError` if no compatible bijection of flags exists.
"""
function build_gkm_connection(data::GKMCombinatorialData{R,V,F}; connection_type = "Algorithmic") where {R,V,F}
  C = F <: AbstractOrbifoldFlagWeight ? QQFieldElem : R
  transport = Dict{Edge, Vector{Int}}()
  a = Dict{Edge, Vector{C}}()

  for unoriented_edge in edges(graph(data))
    e = Edge(src(unoriented_edge), dst(unoriented_edge))
    image, coefficients = _connection_along_edge(data, e)
    inverse_image = invperm(image)
    reverse_coefficients = _coefficients_for_transport(data, reverse(e), inverse_image, C)
    transport[e], transport[reverse(e)] = image, inverse_image
    a[e], a[reverse(e)] = coefficients, reverse_coefficients
  end

  return Connection{C}(transport, a, connection_type)
end

function _coefficients_for_transport(data, e::Edge, image, ::Type{R}) where {R}
  source_edge_flag, _ = _edge_flag_indices(data, e)
  edge_weight = _connection_weight(data, src(e), source_edge_flag)
  result = Vector{R}(undef, length(image))
  for i in eachindex(image)
    difference = _connection_weight(data, src(e), i) -
                 _connection_weight(data, dst(e), image[i])
    coefficient = _connection_coefficient(difference, edge_weight, R)
    isnothing(coefficient) &&
      throw(ArgumentError("transport does not admit a coefficient for flag $i along $e"))
    result[i] = coefficient
  end
  return result
end

# A perfect matching is necessary when a flag has several compatible images.
function _connection_along_edge(data::GKMCombinatorialData{R,V,F}, e::Edge) where {R,V,F}
  C = F <: AbstractOrbifoldFlagWeight ? QQFieldElem : R
  source_flags, target_flags = flags(data, src(e)), flags(data, dst(e))
  length(source_flags) == length(target_flags) ||
    throw(ArgumentError("the endpoints of $e have different numbers of flags"))

  n = length(source_flags)
  candidates = Vector{Vector{Tuple{Int,C}}}(undef, n)
  source_edge_flag, _ = _edge_flag_indices(data, e)
  edge_weight = _connection_weight(data, src(e), source_edge_flag)
  iszero(edge_weight) && throw(ArgumentError("the weight of $e is zero"))

  for i in 1:n
    candidates[i] = Tuple{Int,C}[]
    for j in 1:n
      coefficient = _connection_coefficient(
        _connection_weight(data, src(e), i) - _connection_weight(data, dst(e), j),
        edge_weight,
        C,
      )
      isnothing(coefficient) || push!(candidates[i], (j, coefficient))
    end
    isempty(candidates[i]) &&
      throw(ArgumentError("no compatible image for flag $i along $e"))
  end

  image, coefficients, used = zeros(Int, n), Vector{C}(undef, n), falses(n)
  _find_connection_matching!(image, coefficients, used, candidates, 1) ||
    throw(ArgumentError("the compatible flag pairs along $e do not form a bijection"))
  return image, coefficients
end

_connection_weight(data::GKMCombinatorialData, v::Int, i::Int) = weight(data, v, i)

function _connection_weight(
  data::GKMCombinatorialData{R,V,F}, v::Int, i::Int
) where {R,V,F<:AbstractOrbifoldFlagWeight}
  flag = flags(data, v)[i]
  order = order_of_generic_stabilizer(flag)
  return [QQ(flag.weight[k]) / order for k in 1:rank(lattice(data))]
end

function _find_connection_matching!(image, coefficients, used, candidates, i)
  i > length(candidates) && return true
  for (j, coefficient) in candidates[i]
    used[j] && continue
    image[i], coefficients[i], used[j] = j, coefficient, true
    _find_connection_matching!(image, coefficients, used, candidates, i + 1) && return true
    used[j] = false
  end
  return false
end

function _connection_coefficient(difference, edge_weight, ::Type{R}) where {R}
  indices = edge_weight isa AbstractVector ? eachindex(edge_weight) : 1:rank(parent(edge_weight))
  coordinate = findfirst(i -> !iszero(edge_weight[i]), indices)
  isnothing(coordinate) && return nothing
  coefficient = try
    divexact(difference[coordinate], edge_weight[coordinate])
  catch
    return nothing
  end
  coefficient isa R || return nothing
  difference == coefficient * edge_weight || return nothing
  return coefficient
end

@doc raw"""
    is_valid(G::AbstractGKMGraph, con::AbstractGKMConnection;
             print_diagnostics=false) -> Bool

Check that `con` defines a valid connection on `G`. Every oriented edge must
have a permutation of the flags, inverse transport on the reversed edge, and
coefficients satisfying the connection equation.
"""
function is_valid(G::AbstractGKMGraph, con::AbstractGKMConnection;
                  print_diagnostics::Bool=false)
  trans, coeffs = transport(con), coefficients(con)

  for base_edge in edges(G), e in (base_edge, reverse(base_edge))
    haskey(trans, e) || return _invalid_connection("missing transport for $e", print_diagnostics)
    haskey(coeffs, e) || return _invalid_connection("missing coefficients for $e", print_diagnostics)

    n_source = length(flags(G, src(e)))
    n_target = length(flags(G, dst(e)))
    length(trans[e]) == n_source ||
      return _invalid_connection("transport along $e has the wrong length", print_diagnostics)
    length(coeffs[e]) == n_source ||
      return _invalid_connection("coefficients along $e have the wrong length", print_diagnostics)
    sort(trans[e]) == collect(1:n_target) ||
      return _invalid_connection("transport along $e is not a permutation", print_diagnostics)

    reverse_transport = get(trans, reverse(e), Int[])
    length(reverse_transport) == n_target ||
      return _invalid_connection("reverse transport for $e has the wrong length", print_diagnostics)
    all(reverse_transport[trans[e][i]] == i for i in 1:n_source) ||
      return _invalid_connection("forward and reverse transports for $e are not inverse", print_diagnostics)

    source_edge_flag, target_edge_flag = _edge_flag_indices(G, e)
    trans[e][source_edge_flag] == target_edge_flag ||
      return _invalid_connection("transport along $e does not preserve the edge flag", print_diagnostics)

    for i in 1:n_source
      j, coefficient = trans[e][i], coeffs[e][i]
      data = core(G)
      source_edge_flag, _ = _edge_flag_indices(data, e)
      edge_weight = _connection_weight(data, src(e), source_edge_flag)
      _connection_weight(data, dst(e), j) ==
        _connection_weight(data, src(e), i) - coefficient * edge_weight ||
        return _invalid_connection("connection equation fails along $e at flag $i", print_diagnostics)
    end
  end
  return true
end

is_valid(con::AbstractGKMConnection, G::AbstractGKMGraph; kwargs...) =
  is_valid(G, con; kwargs...)

function _invalid_connection(message::AbstractString, print_diagnostics::Bool)
  print_diagnostics && println(message)
  return false
end
