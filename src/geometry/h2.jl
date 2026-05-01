mutable struct GKM_H2
  edge_lattice::AbstractAlgebra.FPModule{ZZRingElem}
  H2::AbstractAlgebra.FPModule{ZZRingElem}

  edge_to_gen::Dict{Edge,Int}

  quotient::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}

  dual_cone::Cone{QQFieldElem}
  ray_sum::RayVector{QQFieldElem}

  chern::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}
end