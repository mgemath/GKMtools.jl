function weight_class(G, e::Edge)
  return _weight_class(G, e, gens_coeffRing(G))
end

function _weight_class(G, e::Edge, t::Vector{T}) where {T}
  res = zero(t[1])
  w = weight(G, e)
  for i in 1:rank_torus(G)
    res += w[i] * t[i]
  end
  return res
end

function _flag_weight_class(G, v::Int, i::Int, t::Vector{T}) where {T}
  res = zero(t[1])
  w = weight(G, v, i)
  for j in 1:rank_torus(G)
    res += w[j] * t[j]
  end
  return res
end
