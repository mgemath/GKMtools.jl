###############################################################################
#
#   GKMtools
#
###############################################################################

abstract type GKM_graph end
abstract type AbstractGKM_H2 end
abstract type AbstractGKM_connection end
abstract type AbstractGKM_cohomology_ring end

CurveClass_type = AbstractAlgebra.FPModuleElem{ZZRingElem}
GKM_weight_type = Union{ZZRingElem, QQFieldElem}

@attributes mutable struct AbstractGKM_graph{R <: GKM_weight_type}
  g::Graph
  labels::Vector{String}
  weightType::DataType
  M::AbstractAlgebra.Generic.FreeModule{R} # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}} # weight of the T-action
  # This should always be set and can be accessed directly:
  # It should not be changed.
  equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring} # actual type will be Union{Nothing, GKM_cohomology_ring}
  # Use GKM_second_homology() to access this:
  curveClasses::Union{Nothing, AbstractGKM_H2} # actual type will be Union{Nothing, GKM_H2}
  # Use get_GKM_connection() to access this:
  connection::Union{Nothing, AbstractGKM_connection} # actual type will be Union{Nothing, GKM_connection}
  anyConnection::Union{Nothing, AbstractGKM_connection}
  # Use QH_structure_constants() to access this.
  QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}}
  # This should be set to true if the GKM graph is strictly NEF and all relevant
  # QH cohomology structure constants are calculated and stored in QH_structure_consts.
  know_all_QH_structure_consts::Bool
  # Base change matrix for changing to preferred basis for printing quantum products
  QH_preferred_basis::Union{Nothing, AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}}

  function AbstractGKM_graph(
    g::Graph,
    labels::Vector{String},
    M::AbstractAlgebra.Generic.FreeModule{R}, # character group
    w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}},
    equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring},
    curveClasses::Union{Nothing, AbstractGKM_H2},
    connection::Union{Nothing, AbstractGKM_connection},
    QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}},
    know_all_QH_structure_consts::Bool
  ) where R <: GKM_weight_type
    return new{R}(g, labels, R, M, w, equivariantCohomology, curveClasses, connection, nothing, QH_structure_consts, know_all_QH_structure_consts, nothing)
  end
end

@attributes mutable struct AbstractGKM_subgraph
  super::AbstractGKM_graph
  self::AbstractGKM_graph # the GKM subgraph which forgets about the supergraph
  vDict::Vector{Int64} # track how vertices of the subgraph are mapped to that of the supergraph (since Oscar always uses {1, ..., n} as vertex set)

  function AbstractGKM_subgraph(super::AbstractGKM_graph, self::AbstractGKM_graph, vDict::Vector{Int64})
    return new(super, self, vDict)
  end
end

struct GKM_cohomology_ring <: AbstractGKM_cohomology_ring
  gkm::AbstractGKM_graph
  coeffRing::QQMPolyRing # H_T^*(point;Q)
  coeffRingLocalized::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}
  # H_T^*(X; Q), but without checks for consistency (see isGKMclass in cohomology.jl):
  cohomRing::FreeMod{QQMPolyRingElem}
  # H_T^*(X;Q) tensored with the fraction field of H_T^*(point):
  cohomRingLocalized::AbstractAlgebra.Generic.FreeModule{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  edgeWeightClasses::Dict{Edge, QQMPolyRingElem}
  pointEulerClasses::Vector{Union{Nothing, QQMPolyRingElem}}

  function GKM_cohomology_ring(
    gkm::AbstractGKM_graph,
    coeffRing::QQMPolyRing,
    coeffRingLocalized,
    cohomRing::FreeMod{QQMPolyRingElem},
    cohomRingLocalized,
    edgeWeightClasses::Dict{Edge, QQMPolyRingElem},
    pointEulerClasses::Vector{Union{Nothing, QQMPolyRingElem}}
  )
    return new(gkm, coeffRing, coeffRingLocalized, cohomRing, cohomRingLocalized, edgeWeightClasses, pointEulerClasses)
  end
end

mutable struct GKM_connection <: AbstractGKM_connection
  gkm::AbstractGKM_graph
  con::Dict{Tuple{Edge, Edge}, Edge} # assigns to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} # w[e'_i] = w [e_i] - a_i * w[e]

  function GKM_connection(
    gkm::AbstractGKM_graph,
    con::Dict{Tuple{Edge, Edge}, Edge},
    a::Dict{Tuple{Edge, Edge},ZZRingElem}
  )
    return new(gkm, con, a)
  end
end

mutable struct GKM_H2 <: AbstractGKM_H2
  gkm::AbstractGKM_graph
  edgeLattice::AbstractAlgebra.FPModule{ZZRingElem}
  H2::AbstractAlgebra.FPModule{ZZRingElem} # quotient of edgeLattice by relations in H_2.
  edgeToGenIndex::Dict{Edge, Int64}
  quotientMap::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem} # Z-module homomorphism from edgeLattice to H2
  dualConeRaySum::RayVector{QQFieldElem} # sum of rays of the dual cone of the edgeCurveClasses, normalized so that the minimum of pairings with edge curve classes is 1.
  dualCone::Cone{QQFieldElem} # dual cone of the cone of effective curve classes
  chernNumber::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem} # Z-module homomorphism from H2 to ZZ, giving the curve class evaluated on the first chern class of the tangent bundle of the space.
  # The following are used by Seidel spaces to project curve classes to the underlying P1 and the vertical fibres.
  sectionCount::Union{Nothing, AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}}
  verticalProjection::Union{Nothing, AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}}

  function GKM_H2(
    gkm::AbstractGKM_graph,
    edgeLattice::AbstractAlgebra.FPModule{ZZRingElem},
    H2::AbstractAlgebra.FPModule{ZZRingElem},
    edgeToGenIndex::Dict{Edge, Int64},
    quotientMap::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem},
    dualConeRaySum::RayVector{QQFieldElem},
    dualCone::Cone{QQFieldElem},
    chernNumber::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem},
    sectionCount::Union{Nothing, AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}},
    verticalProjection::Union{Nothing, AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}}
  )
    return new(gkm, edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum, dualCone, chernNumber, sectionCount, verticalProjection)
  end
end

@attributes mutable struct GKM_vector_bundle{R <: GKM_weight_type}
  gkm::AbstractGKM_graph{R}
  # Character group of the torus acting on the vector bundle.
  # This could be larger than the torus acting on the gkm space, for example by scaling fibres.
  M::AbstractAlgebra.Generic.FreeModule{R}
  # The homomorphism injecting gkm.M into M.
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}
  # weights[i, j] is the jth weight at the ith vertex.
  w::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}
  # connection along edges of GKM.g for the fibre sub-line-bundles.
  con::Union{Nothing, Dict{Tuple{Edge, Int64}, Int64}}

  function GKM_vector_bundle(
    gkm::AbstractGKM_graph,
    M::AbstractAlgebra.Generic.FreeModule{R},
    GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
    w::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}},
    con::Union{Nothing, Dict{Tuple{Edge, Int64}, Int64}}
  ) where R <: GKM_weight_type
    return new{R}(gkm, M, GMtoM, w, con)
  end
end