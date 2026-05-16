struct Vertex <: AbstractVertex
  label::String
end

struct ToricVertex <: AbstractVertex
  label::String
end

struct InertiaVertex{T} <: AbstractVertex where {T <: AbstractVertex}
  label::String
  old_vertex::T
end

function get_string(V::AbstractVertex)
  return V.label
end