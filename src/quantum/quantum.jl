mutable struct GKM_Quantum
  structure_constants::Dict{Any,Array{Any,3}}

  known_complete::Bool

  preferred_basis::Union{Nothing,Any}
end