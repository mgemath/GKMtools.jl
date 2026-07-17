struct OrbifoldGKMVectorBundle{
  R,
  V,
  F,
  X<:AbstractOrbifoldGKMGraph{R,V,F},
} <: AbstractOrbifoldGKMVectorBundle{R,V,F}
  base::X
  M::AbstractAlgebra.Generic.FreeModule{R}
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}
  weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}
  fiber_reps::Vector{ZZMatrix}
  connection::Connection{R}
end

function Base.show(io::IO, E::OrbifoldGKMVectorBundle)
  print(
    io,
    "Orbifold GKM vector bundle of rank $(rank(E)) over $(num_vertices(baseof(E))) stacky vertices",
  )
  show(io, MIME"text/plain"(), E)
end

function Base.show(io::IO, ::MIME"text/plain", E::OrbifoldGKMVectorBundle)
  println(
    io,
    "Orbifold GKM vector bundle of rank $(rank(E)) over $(num_vertices(baseof(E))) stacky vertices with weights:",
  )
  for v in 1:num_vertices(baseof(E))
    print(io, "$(label(baseof(E), v)): ")
    print(io, join([string(E.weights[v, i]) for i in 1:rank(E)], ", "))
    if !isempty(baseof(E).vertex_isotropy[v].isotropy_group)
      print(io, " | fibre representation ")
      show(io, MIME"text/plain"(), E.fiber_reps[v])
    end
    v == num_vertices(baseof(E)) || println(io)
  end
end

struct GKMVectorBundle{
  R,
  V,
  F,
  X<:AbstractGKMGraph{R,V,F},
} <: AbstractGKMVectorBundle{R,V,F}
  base::X
  M::AbstractAlgebra.Generic.FreeModule{R}
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}
  weights::Matrix{
    AbstractAlgebra.Generic.FreeModuleElem{R}
  }
  connection::Connection{R}
end

function Base.show(io::IO, E::GKMVectorBundle)
  print(
    io,
    "GKM vector bundle of rank $(rank(E)) over " *
    "$(num_vertices(baseof(E))) vertices",
  )
end

function Base.show(
  io::IO,
  ::MIME"text/plain",
  E::GKMVectorBundle,
)
  println(
    io,
    "GKM vector bundle of rank $(rank(E)) over " *
    "$(num_vertices(baseof(E))) vertices with weights:",
  )

  for v in 1:num_vertices(baseof(E))
    print(io, "$(label(baseof(E), v)): ")
    print(
      io,
      join(
        (string(E.weights[v, i]) for i in 1:rank(E)),
        ", ",
      ),
    )

    v == num_vertices(baseof(E)) || println(io)
  end
end
function _edge_signature(G::AbstractGKMGraph)
  return sort([
    (min(src(e), dst(e)), max(src(e), dst(e)))
    for e in edges(G)
  ])
end

function _same_flag_data(G::AbstractGKMGraph, H::AbstractGKMGraph, v::Int)
  length(flags(G, v)) == length(flags(H, v)) || return false

  for i in eachindex(flags(G, v))
    weight(G, v, i) == weight(H, v, i) || return false

    flag_G = flags(G, v)[i]
    flag_H = flags(H, v)[i]
    if flag_G isa AbstractOrbifoldFlagWeight || flag_H isa AbstractOrbifoldFlagWeight
      flag_G isa AbstractOrbifoldFlagWeight || return false
      flag_H isa AbstractOrbifoldFlagWeight || return false
      order_of_generic_stabilizer(flag_G) == order_of_generic_stabilizer(flag_H) || return false
    end
  end

  return true
end

function _same_graph_data(G::AbstractGKMGraph, H::AbstractGKMGraph)
  num_vertices(G) == num_vertices(H) || return false
  rank_torus(G) == rank_torus(H) || return false
  [label(G, v) for v in vertices(G)] == [label(H, v) for v in vertices(H)] || return false
  _edge_signature(G) == _edge_signature(H) || return false

  for v in vertices(G)
    _same_flag_data(G, H, v) || return false
  end

  return true
end

function _same_orbifold_isotropy_data(G::AbstractOrbifoldGKMGraph, H::AbstractOrbifoldGKMGraph)
  length(G.vertex_isotropy) == length(H.vertex_isotropy) || return false
  length(G.flag_isotropy) == length(H.flag_isotropy) || return false

  for v in vertices(G)
    G.vertex_isotropy[v].isotropy_group == H.vertex_isotropy[v].isotropy_group || return false
    G.vertex_isotropy[v].tangent_rep == H.vertex_isotropy[v].tangent_rep || return false
    length(G.flag_isotropy[v]) == length(H.flag_isotropy[v]) || return false

    for i in eachindex(G.flag_isotropy[v])
      G.flag_isotropy[v][i].isotropy_group == H.flag_isotropy[v][i].isotropy_group || return false
      G.flag_isotropy[v][i].embedding == H.flag_isotropy[v][i].embedding || return false
    end
  end

  return true
end

function Base.:(==)(G::GKMGraph, H::GKMGraph)
  return _same_graph_data(G, H)
end

function Base.:(==)(G::OrbifoldGKMGraph, H::OrbifoldGKMGraph)
  return (
    _same_graph_data(G, H) &&
    _same_orbifold_isotropy_data(G, H)
  )
end

function Base.:(==)(E::GKMVectorBundle, F::GKMVectorBundle)
  return (
    _same_graph_data(baseof(E), baseof(F)) &&
    E.M == F.M &&
    E.GMtoM == F.GMtoM &&
    E.weights == F.weights
  )
end

function Base.:(==)(E::OrbifoldGKMVectorBundle, F::OrbifoldGKMVectorBundle)
  return (
    _same_graph_data(baseof(E), baseof(F)) &&
    _same_orbifold_isotropy_data(baseof(E), baseof(F)) &&
    E.M == F.M &&
    E.GMtoM == F.GMtoM &&
    E.weights == F.weights &&
    E.fiber_reps == F.fiber_reps
  )
end
