struct GKMFlag{R}
  vertex::Int
  weight::AbstractAlgebra.Generic.FreeModuleElem{R}
  edge::Union{Nothing,Edge}
end