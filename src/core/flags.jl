struct FlagWeight{R} <: AbstractFlagWeight{R}
  weight::AbstractAlgebra.Generic.FreeModuleElem{R}
end

struct ToricFlagWeight{R} <: AbstractFlagWeight{R}
  weight::AbstractAlgebra.Generic.FreeModuleElem{R}
end

struct OrbifoldFlagWeight{R} <: AbstractOrbifoldFlagWeight{R}
  weight::AbstractAlgebra.Generic.FreeModuleElem{R} # This is the primitive weight of the flag, which may be an integer multiple of the actual weight if we are in the orbifold setting
  order_of_generic_stabilizer::Int # This is the order of the generic stabilizer along the edge corresponding to this flag, which we can use to recover the actual weight in the orbifold setting
end

struct OrbifoldToricFlagWeight{R} <: AbstractOrbifoldFlagWeight{R}
  weight::AbstractAlgebra.Generic.FreeModuleElem{R} # This is the primitive weight of the flag, which may be an integer multiple of the actual weight if we are in the orbifold setting
  order_of_generic_stabilizer::Int # This is the order of the generic stabilizer along the edge corresponding to this flag, which we can use to recover the actual weight in the orbifold setting
end

order_of_generic_stabilizer(fw::AbstractOrbifoldFlagWeight) = fw.order_of_generic_stabilizer
