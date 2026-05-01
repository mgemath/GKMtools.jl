mutable struct GKM_Cohomology
  coeff_ring::QQMPolyRing
  localized_ring::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}

  cohomology::FreeMod{QQMPolyRingElem}

  localized_cohomology::AbstractAlgebra.Generic.FreeModule{
    AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}
  }

  edge_classes::Dict{Edge,QQMPolyRingElem}

  euler_classes::Vector{Union{Nothing,QQMPolyRingElem}}
end