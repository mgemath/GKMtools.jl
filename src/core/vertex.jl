struct Vertex <: AbstractVertex
  label::String
end

## Toric specific types

struct ToricVertex <: AbstractVertex
  label::String
end

## Flag specific types
abstract type AbstractFlagVertex <: AbstractVertex end
struct FlagVertex{L} <: AbstractFlagVertex where {L}
  label::String
  flag::L
end

struct GeneralizedFlagVertex <: AbstractFlagVertex
  label::String
  flag::WeylGroupElem
end

flag(V::AbstractFlagVertex) = V.flag
# Base.length(V::GeneralizedFlagVertex) = length(flag(V))
Base.isless(V::GeneralizedFlagVertex, W::GeneralizedFlagVertex) = flag(V) < flag(W)
Base.isless(V::FlagVertex, W::FlagVertex) = error("Not implemented")
struct BlowupVertex{T} <: AbstractVertex where {T <: AbstractVertex}
  label::String
  old_vertex::T
end

struct ProductVertex{V1<:AbstractVertex,V2<:AbstractVertex} <: AbstractVertex
  first::V1
  second::V2
end

first_vertex(V::ProductVertex) = V.first
second_vertex(V::ProductVertex) = V.second
get_string(v::ProductVertex) = "$(get_string(v.first)),$(get_string(v.second))"

struct InertiaVertex{T} <: AbstractVertex where {T <: AbstractVertex}
  label::String
  old_vertex::T
end

function get_string(V::AbstractVertex)
  return V.label
end