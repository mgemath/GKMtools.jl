###############################################################################
#
#   GKMgraphs
#
###############################################################################

abstract type AbstractGKM_graph end
abstract type AbstractGKM_graph_w_connection <: AbstractGKM_graph end
abstract type AbstractGKM_subgraph end
# abstract type AbstractGKM_H2 end
# abstract type AbstractGKM_connection end
# abstract type AbstractGKM_cohomology_ring end

# CurveClass_type = AbstractAlgebra.FPModuleElem{ZZRingElem}
# GKM_weight_type = Union{ZZRingElem, QQFieldElem}
GKM_weight_type = Union{Int64, Bool}

@attributes mutable struct GKM_graph{R <: GKM_weight_type} <: AbstractGKM_graph
  g::Graph
  labels::Vector{String}
  # M::AbstractAlgebra.Generic.FreeModule{R} # character group
  # w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}} # weight of the T-action
  # This should always be set and can be accessed directly:
  # It should not be changed.
  # equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring} # actual type will be Union{Nothing, GKM_cohomology_ring}
  # Use GKM_second_homology() to access this:
  # curveClasses::Union{Nothing, AbstractGKM_H2} # actual type will be Union{Nothing, GKM_H2}
  # Use get_GKM_connection() to access this:
  connection::Union{Nothing, AbstractGKM_connection} # actual type will be Union{Nothing, GKM_connection}
  # Use QH_structure_constants() to access this.
  # QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}}
  # This should be set to true if the GKM graph is strictly NEF and all relevant
  # QH cohomology structure constants are calculated and stored in QH_structure_consts.
  know_all_QH_structure_consts::Bool

  function GKM_graph(
    g::Graph,
    labels::Vector{String},
    # M::AbstractAlgebra.Generic.FreeModule{R}, # character group
    # w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}},
    # equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring},
    # curveClasses::Union{Nothing, AbstractGKM_H2},
    # connection::Union{Nothing, AbstractGKM_connection},
    # QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}},
    know_all_QH_structure_consts::Bool
  ) where R <: GKM_weight_type
    return new{R}(g, labels, know_all_QH_structure_consts)
    # return new{R}(g, labels, R, M, w, equivariantCohomology, curveClasses, connection, QH_structure_consts, know_all_QH_structure_consts)
  end
end

mutable struct GKM_graph_w_connection{R <: GKM_weight_type} <: AbstractGKM_graph_w_connection
  g::Graph
  labels::Vector{String}
  # M::AbstractAlgebra.Generic.FreeModule{R} # character group
  # w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}} # weight of the T-action
  # This should always be set and can be accessed directly:
  # It should not be changed.
  # equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring} # actual type will be Union{Nothing, GKM_cohomology_ring}
  # Use GKM_second_homology() to access this:
  # curveClasses::Union{Nothing, AbstractGKM_H2} # actual type will be Union{Nothing, GKM_H2}
  con::Dict{Tuple{Edge, Edge}, Edge} # assigns to each edges e & e_i with src(e)=src(e_i) an edge e'_i with src(e'_i)=dst(e).
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} # w[e'_i] = w [e_i] - a_i * w[e]
  # Use QH_structure_constants() to access this.
  # QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}}
  # This should be set to true if the GKM graph is strictly NEF and all relevant
  # QH cohomology structure constants are calculated and stored in QH_structure_consts.
  know_all_QH_structure_consts::Bool

  function GKM_graph(
    g::Graph,
    labels::Vector{String},
    # M::AbstractAlgebra.Generic.FreeModule{R}, # character group
    # w::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{R}},
    # equivariantCohomology::Union{Nothing, AbstractGKM_cohomology_ring},
    # curveClasses::Union{Nothing, AbstractGKM_H2},
    # connection::Union{Nothing, AbstractGKM_connection},
    # QH_structure_consts::Dict{CurveClass_type, Array{Any, 3}},
    know_all_QH_structure_consts::Bool
  ) where R <: GKM_weight_type
    return new{R}(g, labels, know_all_QH_structure_consts)
    # return new{R}(g, labels, R, M, w, equivariantCohomology, curveClasses, connection, QH_structure_consts, know_all_QH_structure_consts)
  end
end