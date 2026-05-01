struct ConnectionData
  transport::Dict{Edge,Vector{Int}}
  a::Dict{Edge,Vector{ZZRingElem}}
end

struct CanonicalConnection <: GKMConnection
  data::ConnectionData
end

struct AlgorithmicConnection <: GKMConnection
  data::ConnectionData
end

# Interface

transport(c::GKMConnection) = c.data.transport
coefficients(c::GKMConnection) = c.data.a

is_canonical(::CanonicalConnection) = true
is_canonical(::GKMConnection) = false