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