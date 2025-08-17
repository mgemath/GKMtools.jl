###############################################################################
#
#   GKMtools
#
###############################################################################

struct GKM_cohomology_ring <: AbstractGKM_cohomology_ring
  coeffRing::QQMPolyRing # H_T^*(point;Q)
  coeffRingLocalized::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}
  # H_T^*(X; Q), but without checks for consistency (see isGKMclass in cohomology.jl):
  cohomRing::FreeMod{QQMPolyRingElem}
  # H_T^*(X;Q) tensored with the fraction field of H_T^*(point):
  cohomRingLocalized::AbstractAlgebra.Generic.FreeModule{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
  edgeWeightClasses::Dict{Edge, QQMPolyRingElem}
  pointEulerClasses::Vector{Union{Nothing, QQMPolyRingElem}}

  function GKM_cohomology_ring(
    coeffRing::QQMPolyRing,
    coeffRingLocalized,
    cohomRing::FreeMod{QQMPolyRingElem},
    cohomRingLocalized,
    edgeWeightClasses::Dict{Edge, QQMPolyRingElem},
    pointEulerClasses::Vector{Union{Nothing, QQMPolyRingElem}}
  )
    return new(coeffRing, coeffRingLocalized, cohomRing, cohomRingLocalized, edgeWeightClasses, pointEulerClasses)
  end
end

mutable struct GKM_connection <: AbstractGKM_connection
  con::Dict{Tuple{Edge, Edge}, Edge} # assigns to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} # w[e'_i] = w [e_i] - a_i * w[e]

  function GKM_connection(
    con::Dict{Tuple{Edge, Edge}, Edge},
    a::Dict{Tuple{Edge, Edge}, ZZRingElem}
  )
    return new(con, a)
  end
end

mutable struct GKM_H2 <: AbstractGKM_H2
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
    return new(edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum, dualCone, chernNumber, sectionCount, verticalProjection)
  end
end

@attributes mutable struct GKM_graph{R <: Weight} <: AbstractGKM_graph{R}
  g::Graph
  labels::Vector{String}
  weightType::DataType
  M::AbstractAlgebra.Generic.FreeModule{R} # character group
  w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}} # weight of the T-action
  # This should always be set and can be accessed directly:
  # It should not be changed.
  equivariantCohomology::Union{Nothing, GKM_cohomology_ring} # actual type will be Union{Nothing, GKM_cohomology_ring}
  # # Use GKM_second_homology() to access this:
  H2::Union{Nothing, GKM_H2} # actual type will be Union{Nothing, GKM_H2}
  # # Use get_GKM_connection() to access this:
  connection::Union{Nothing, GKM_connection} # actual type will be Union{Nothing, GKM_connection}
  # Use QH_structure_constants() to access this.
  QH_structure_consts::Dict{CurveClass, Array{Any, 3}}
  # This should be set to true if the GKM graph is strictly NEF and all relevant
  # QH cohomology structure constants are calculated and stored in QH_structure_consts.
  know_all_QH_structure_consts::Bool
  # Base change matrix for changing to preferred basis for printing quantum products
  QH_preferred_basis::Union{Nothing, AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}}

  function GKM_graph(
    g::Graph,
    labels::Vector{String},
    M::AbstractAlgebra.Generic.FreeModule{R}, # character group
    w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}},
    equivariantCohomology::Union{Nothing, GKM_cohomology_ring},
    H2::Union{Nothing, GKM_H2},
    connection::Union{Nothing, GKM_connection},
    QH_structure_consts::Dict{CurveClass, Array{Any, 3}},
    know_all_QH_structure_consts::Bool
  ) where R <: Weight
    return new{R}(g, labels, R, M, w, equivariantCohomology, H2, connection, QH_structure_consts, know_all_QH_structure_consts, nothing)
    # return new{R}(g, labels, R, M, w, equivariantCohomology, H2, connection, QH_structure_consts, know_all_QH_structure_consts, nothing)
  end
end

mutable struct GKM_morphism{R <: Weight} <: AbstractGKM_morphism{R}
  Domain::GKM_graph{R}
  Codomain::GKM_graph{R} # the GKM subgraph which forgets about the supergraph
  vDict::Vector{Int64} # track how vertices of the subgraph are mapped to that of the supergraph (since Oscar always uses {1, ..., n} as vertex set)

  function GKM_morphism(Domain::GKM_graph{R}, Codomain::GKM_graph{R}, vDict::Vector{Int64}) where R <: Weight
    return new{R}(Domain, Codomain, vDict)
  end
end

mutable struct GKM_subgraph{R <: Weight} <: AbstractGKM_subgraph{R}   # an injective morphism og GKM graphs
  Domain::GKM_graph{R}
  Codomain::GKM_graph{R} # the GKM subgraph which forgets about the supergraph
  vDict::Vector{Int64} # track how vertices of the subgraph are mapped to that of the supergraph (since Oscar always uses {1, ..., n} as vertex set)

  function GKM_subgraph(Domain::GKM_graph{R}, Codomain::GKM_graph{R}, vDict::Vector{Int64}) where R <: Weight
    return new{R}(Domain, Codomain, vDict)
  end
end

# @attributes mutable struct GKM_vector_bundle{R <: Weight}
#   gkm::AbstractGKM_graph{R}
#   # Character group of the torus acting on the vector bundle.
#   # This could be larger than the torus acting on the gkm space, for example by scaling fibres.
#   M::AbstractAlgebra.Generic.FreeModule{R}
#   # The homomorphism injecting gkm.M into M.
#   GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}
#   # weights[i, j] is the jth weight at the ith vertex.
#   w::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}
#   # connection along edges of GKM.g for the fibre sub-line-bundles.
#   con::Union{Nothing, Dict{Tuple{Edge, Int64}, Int64}}

#   function GKM_vector_bundle(
#     gkm::AbstractGKM_graph,
#     M::AbstractAlgebra.Generic.FreeModule{R},
#     GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
#     w::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}},
#     con::Union{Nothing, Dict{Tuple{Edge, Int64}, Int64}}
#   ) where R <: Weight
#     return new{R}(gkm, M, GMtoM, w, con)
#   end
# end