struct Vertex <: AbstractVertex
  label::String
end

struct ToricVertex <: AbstractVertex
  label::String
end

struct FlagVertex{L} <: AbstractVertex where {L}
  label::String
  flag::L
end

struct GeneralizedFlagVertex <: AbstractVertex
  label::String
  flag::WeylGroupElem
end

struct BlowupVertex{T} <: AbstractVertex where {T <: AbstractVertex}
  label::String
  old_vertex::T
end

struct InertiaVertex{T} <: AbstractVertex where {T <: AbstractVertex}
  label::String
  old_vertex::T
end

function get_string(V::AbstractVertex)
  return V.label
end