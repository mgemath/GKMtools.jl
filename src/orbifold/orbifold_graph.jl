# This file is part of GKMtools.jl, licensed under the MIT License (MIT).
struct OrbifoldGKMGraph{R,V,F} <: AbstractOrbifoldGKMGraph{R,V,F}
  core::GKMCombinatorialData{R,V,F}
  vertex_isotropy::Vector{OrbifoldVertexIsotropy} # isotropy data for each vertex
  flag_isotropy::Vector{Vector{OrbifoldFlagIsotropy}} # isotropy data for each flag
end

# Interface delegation

# function weight(g::OrbifoldGKMGraph, e::Edge)
#   return g.edge_multiplicity[e] * weight(g.core, e)
# end


# --- Isotropy Structures ---

function Base.show(io::IO, I::OrbifoldVertexIsotropy)
  print(io, "VertexIsotropy(order $(order_of_isotropy_group(I)))")
end

function Base.show(io::IO, ::MIME"text/plain", I::OrbifoldVertexIsotropy)
  println(io, "Orbifold Vertex Isotropy:")
  println(io, "  Group Structure (cyclic factors): ", I.isotropy_group)
  println(io, "  Tangent Representation: ")
  show(io, "text/plain", I.tangent_rep)
end

function Base.show(io::IO, I::OrbifoldFlagIsotropy)
  print(io, "FlagIsotropy(order $(order_of_isotropy_group(I)))")
end

function Base.show(io::IO, ::MIME"text/plain", I::OrbifoldFlagIsotropy)
  println(io, "Orbifold Flag Isotropy:")
  println(io, "  Group Structure: ", I.isotropy_group)
  println(io, "  Embedding Matrix: ")
  show(io, "text/plain", I.embedding)
end

# --- OrbifoldGKMGraph ---

function Base.show(io::IO, G::OrbifoldGKMGraph)
  print(io, "OrbifoldGKMGraph with $(nv(G.core.g)) stacky vertices")
end

function Base.show(io::IO, ::MIME"text/plain", G::OrbifoldGKMGraph)
  println(io, "Orbifold GKM Graph")
  println(io, "---------------------------")
  show(io, MIME"text/plain"(), G.core)

  # Summary of stacky data
  v_orders = [order_of_isotropy_group(v) for v in G.vertex_isotropy]
  println(io, "  Vertex Isotropy Orders: ", v_orders)

  # Example of printing a specific vertex detail if graph is small
  if nv(G.core.g) <= 5
    for i in 1:nv(G.core.g)
      println(
        io,
        "  Vertex $i: Order $(v_orders[i]), $(length(G.flag_isotropy[i])) incident flags",
      )
    end
  end
end