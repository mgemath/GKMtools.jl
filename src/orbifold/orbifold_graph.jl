# This file is part of GKMtools.jl, licensed under the MIT License (MIT).
struct OrbifoldGKMGraph{R,V,F} <: AbstractOrbifoldGKMGraph{R,V,F}
  core::GKMCombinatorialData{R,V,F}
  vertex_isotropy::Vector{OrbifoldVertexIsotropy} # isotropy data for each vertex
  flag_isotropy::Vector{Vector{OrbifoldFlagIsotropy}} # isotropy data for each flag

  cohomology::GKMCohomology # cohomology of the orbifold GKM graph

  connection::Connection
  # H2::GKM_H2
  H2::Nothing ## TODO: implement H2 for orbifolds



  function OrbifoldGKMGraph(core::GKMCombinatorialData{R,V,F}, vertex_isotropy::Vector{OrbifoldVertexIsotropy}, flag_isotropy::Vector{Vector{OrbifoldFlagIsotropy}}, connection::Connection) where {R,V,F}
    # Check that the number of vertices matches
    nv(core.g) == length(vertex_isotropy) || throw(ArgumentError("Number of vertices in core graph does not match length of vertex_isotropy"))
    nv(core.g) == length(flag_isotropy) || throw(ArgumentError("Number of vertices in core graph does not match length of flag_isotropy"))
    # for v in 1:nv(core.g)
    #     degree(core.g, v) == length(flag_isotropy[v]) || throw(ArgumentError("Degree of vertex $v does not match length of flag_isotropy[$v]"))
    # end

    # Create the cohomology ring for the orbifold GKM graph
    cohomology = create_cohomology(rank(core.M), nv(core.g))
    t = gens_coeffRing(cohomology.localized_cohomology)
    cohomology.euler_classes = [_euler_class(core, i, t) for i in 1:nv(core.g)]
    # H2 = _GKM_second_homology(core)
    H2 = nothing
    return new{R,V,F}(core, vertex_isotropy, flag_isotropy, cohomology, connection, H2)
  end
end

# Interface delegation

# function weight(g::OrbifoldGKMGraph, e::Edge)
#   return g.edge_multiplicity[e] * weight(g.core, e)
# end

# --- OrbifoldGKMGraph ---

function Base.show(io::IO, G::OrbifoldGKMGraph)
  print(io, "OrbifoldGKMGraph with $(nv(G.core.g)) stacky vertices")
end

function Base.show(io::IO, ::MIME"text/plain", G::OrbifoldGKMGraph)
  print(
    io, "Orbifold GKM graph with $(n_vertices(graph(G))) nodes, valency $(valency(G)) and axial function:"
  )
  for e in edges(G)
    print(io, "\n$(label(G, src(e))) -> $(label(G, dst(e))) => $(weight(G, e))")
  end
  
  if !is_compact(G)
    print(io, "\nStandalone flags:")
  end

  print(io, "\nVertex Isotropy:")
  for v in 1:n_vertices(G.core.g)
    issmooth(G.vertex_isotropy[v]) && continue
    print(io, "\n$(label(G, v)) => ")
    show(io, MIME"text/plain"(), G.vertex_isotropy[v])
  end
  print(io, "\nFlag Isotropy:")
  for v in 1:n_vertices(G.core.g)
    for (i, isotropy) in enumerate(G.flag_isotropy[v])
      issmooth(isotropy) && continue
      print(io, "\n$(label(G, v)).$i => ")
      show(io, MIME"text/plain"(), isotropy)
    end
  end

  print_connection(io, G; verbose = false)
end
#   println(io, "Orbifold GKM Graph")
#   println(io, "---------------------------")
#   show(io, MIME"text/plain"(), G.core)

#   # Summary of stacky data
#   v_orders = [order_of_isotropy_group(v) for v in G.vertex_isotropy]
#   println(io, "  Vertex Isotropy Orders: ", v_orders)

#   # Example of printing a specific vertex detail if graph is small
#   if nv(G.core.g) <= 5
#     for i in 1:nv(G.core.g)
#       println(
#         io,
#         "  Vertex $i: Order $(v_orders[i]), $(length(G.flag_isotropy[i])) incident flags",
#       )
#     end
#   end
# end

@doc raw"""
    convert(OrbifoldGKMGraph{R,V,F}, G::GKMGraph{R,V,F})

Regard a smooth GKM graph as an orbifold GKM graph with trivial isotropy.
The combinatorial data and connection of `G` are reused.
"""
function Base.convert(
  ::Type{OrbifoldGKMGraph{R,V,F}},
  G::GKMGraph{R,V,F},
) where {R,V,F}
  vertex_isotropy = [
    smooth_orbifold_vertex_isotropy_group(length(flags(G, v)))
    for v in vertices(G)
  ]
  flag_isotropy = [
    [
      smooth_orbifold_flag_isotropy_group(length(flags(G, v)), 0)
      for _ in flags(G, v)
    ]
    for v in vertices(G)
  ]

  return OrbifoldGKMGraph(G.core, vertex_isotropy, flag_isotropy, G.connection)
end