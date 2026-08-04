mutable struct GKMCohomology <: AbstractCohomology
  coefficient_ring::QQMPolyRing
  localized_coefficient_ring::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}
  n_vertices::Int
  edge_classes::Dict{Edge,QQMPolyRingElem}
  euler_classes::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
end

function Base.show(io::IO, H::GKMCohomology)
  println(io, "GKM cohomology ring with $(H.n_vertices) vertices")
  print(io, "Coefficient ring: $(H.coefficient_ring)")
  # print(io, "Localized coefficient ring: $(H.localized_coefficient_ring)")
end

"""A GKM cohomology class, stored by its fixed-point restrictions."""
struct GKMClass{G,R}
  parent::GKMCohomology
  graph::G
  restrictions::Vector{R}
end

graph(c::GKMClass) = c.graph


function Base.show(io::IO, c::GKMClass)
  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM class")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM class with $(length(c.restrictions)) nodes")
  end
end

function Base.show(io::IO, ::MIME"text/plain", c::GKMClass)
  print(io, "GKM class with restrictions: \n[", c.restrictions[1])
  for r in c.restrictions[2:end]
    print(io, ", ", r)
  end
  print(io, "]")
end

function get_cohomology(G::AbstractGKMGraph)
  return G.cohomology
end

AbstractAlgebra.coefficient_ring(H::GKMCohomology) = H.coefficient_ring


function create_cohomology(rk_torus, n_vertices_G)
  H_pt, _ = polynomial_ring(QQ, vcat(["t$i" for i in 1:rk_torus]))
  H_pt_loc = fraction_field(H_pt)

  return GKMCohomology(
    H_pt,
    H_pt_loc,
    n_vertices_G,
    Dict{Edge,QQMPolyRingElem}(),
    QQMPolyRingElem[],
  )
end

