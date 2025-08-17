###############################################################################
#
#   GKMtools
#
###############################################################################

Weight = Union{ZZRingElem, QQFieldElem}  # weights are element of Z^n or Q^n
CurveClass = AbstractAlgebra.FPModuleElem{ZZRingElem}

abstract type AbstractGKM_graph{R <: Weight} end  # in the abstract GKM graph, we chose the type of weights
abstract type AbstractGKM_connection end
abstract type AbstractGKM_cohomology_ring end

abstract type AbstractGKM_morphism{R <: Weight} end
abstract type AbstractGKM_subgraph{R <: Weight} <: AbstractGKM_morphism{R} end

abstract type AbstractGKM_H2 end
