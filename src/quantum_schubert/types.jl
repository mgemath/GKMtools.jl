const _QSDegree = Tuple{Vararg{Int}}
const _QSKey = Tuple{Int,Int,Int,_QSDegree}

struct _QSEdge
  vertex::Int
  degree::_QSDegree
  coefficients::Vector{Int}
end

"""
Root data and a lazy cache for ordinary small quantum Schubert calculus.
Construct with [`quantum_schubert_context`](@ref). Basis indices are vertex
indices; Novikov coordinates correspond to `omitted_roots`. Mutable caches
are not thread-safe. Auxiliary equivariant values use QQ or GF(p), selected by the constructor.
"""
struct QuantumSchubertContext{K,T}
  representatives::Vector
  codimensions::Vector{Int}
  omitted_roots::Vector{Int}
  novikov_degrees::Vector{Int}
  unit::Int
  outgoing::Vector{Vector{_QSEdge}}
  incoming::Vector{Vector{_QSEdge}}
  field::K
  prime::Union{Nothing,BigInt}
  diagonal::Vector{T}
  self_restrictions::Vector{T}
  inversion_values::Vector{Vector{T}}
  restriction_cache::Dict{Tuple{Int,Int},T}
  cache::Dict{_QSKey,T}
  diagonal_cache::Dict{Tuple{Int,_QSDegree},T}
  max_cache_entries::Int
  lower_intervals::Union{Nothing,Vector{BitVector}}
  degree_vertices::Vector{Vector{Int}}
end

_qs_zero_degree(Q) = Tuple(zeros(Int, length(Q.omitted_roots)))
_qs_degree_value(Q, d) = sum((a*b for (a,b) in zip(Q.novikov_degrees, d)); init=0)

function _qs_degree(Q, degree)
  d = degree isa Integer ? (degree == 0 ? _qs_zero_degree(Q) : (Int(degree),)) : Tuple(Int.(degree))
  length(d) == length(Q.omitted_roots) || throw(ArgumentError("wrong number of Novikov coordinates"))
  all(>=(0), d) || throw(ArgumentError("Novikov coordinates must be nonnegative"))
  return d
end

function _qs_vertex(Q, u::Integer)
  1 <= u <= length(Q.representatives) || throw(BoundsError(Q.representatives, u))
  return Int(u)
end

function _qs_vertex(Q, u::AbstractString)
  i = findfirst(w -> string(w) == u, Q.representatives)
  isnothing(i) && throw(ArgumentError("unknown Schubert representative: $u"))
  return i
end

function _qs_degrees(Q, bound::Int)
  result = _QSDegree[]
  function visit(prefix, i, remaining)
    if i > length(Q.novikov_degrees)
      push!(result, Tuple(prefix))
      return
    end
    for d in 0:div(remaining, Q.novikov_degrees[i])
      visit(vcat(prefix, d), i+1, remaining-d*Q.novikov_degrees[i])
    end
  end
  bound >= 0 && visit(Int[], 1, bound)
  return result
end

# A generational cache bounds retained entries, including during one product.
# Eviction only discards memoized answers; it cannot change the result.
function _qs_remember!(cache, key, value, limit)
  limit == 0 && return value
  length(cache) >= limit && empty!(cache)
  cache[key] = value
  return value
end

function _qs_field(p)
  isnothing(p) && return QQ
  p isa Integer || throw(ArgumentError("p must be a prime integer"))
  p > 1 && is_prime(ZZ(p)) || throw(ArgumentError("p must be prime"))
  return p <= typemax(Int) ? Oscar.Nemo.Native.GF(Int(p)) : Oscar.Nemo.Native.GF(ZZ(p))
end

_qs_integer(Q, c) = isnothing(Q.prime) ? BigInt(ZZ(c)) : BigInt(lift(ZZ,c))
