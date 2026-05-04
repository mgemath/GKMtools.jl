struct GKMFlag{R}
  vertex::Int
  weight::AbstractAlgebra.Generic.FreeModuleElem{R}
  edge::Union{Nothing,Edge}
end

struct Flag
  vertex::Int
  edge::Union{Nothing,Edge}
end