struct ConnectionData{R}
  transport::Dict{Edge, Vector{Int}}
  a::Dict{Edge, Vector{R}}
end

struct CanonicalConnection{R} <: AbstractGKMConnection{R}
  data::ConnectionData{R}
end

struct AlgorithmicConnection{R} <: AbstractGKMConnection{R}
  data::ConnectionData{R}
end

# Interface

transport(c::AbstractGKMConnection{R}) = c.data.transport
coefficients(c::AbstractGKMConnection{R}) = c.data.a

is_canonical(::CanonicalConnection{R}) = true
is_canonical(::AbstractGKMConnection{R}) = false