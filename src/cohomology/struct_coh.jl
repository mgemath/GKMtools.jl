mutable struct GKMCohomology <: AbstractCohomology
  coefficient_ring::QQMPolyRing
  localized_coefficient_ring::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}
  # coeff_ring::QQMPolyRing
  # localized_ring::AbstractAlgebra.Generic.FracField{QQMPolyRingElem}

  # cohomology::FreeMod{QQMPolyRingElem}
  # cohomology::MPolyQuoRing{QQMPolyRingElem}

  localized_cohomology::MPolyQuoRing{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}}

  edge_classes::Dict{Edge,QQMPolyRingElem}

  euler_classes::Vector{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}
end

function get_cohomology(G::AbstractGKMGraph)
  return G.cohomology
end

function gens_coeffRing(R::MPolyQuoRing{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}})
  return gens(coefficient_ring(R))
end

function gens_coeffRing(G)
  R = get_cohomology(G)
  return gens(coefficient_ring(R.localized_cohomology))
end

function gens_cohomRing(R::GKMCohomology)
  return gens(R.localized_cohomology)
end

function gens_cohomRing(G::AbstractGKMGraph)
  R = get_cohomology(G)
  return gens(R.localized_cohomology)
end

function create_relations(t)
  len_t = length(t)
  len = 1 + div(len_t * (len_t - 1), 2)

  ans = similar(t, len)

  for i in 1:len_t
    for j in i+1:len_t
      ans[div((i - 1) * (2 * len_t - i), 2) + j - i] = t[i] * t[j]
    end
  end

  ans[end] = sum(t) - one(t[1])

  return ans
end

function create_cohomology(rk_torus, n_vertices_G)

  # creation of the localized cohomology ring
  H_pt, _ = polynomial_ring(QQ, vcat(["t$i" for i in 1:rk_torus]))
  H_pt_loc = fraction_field(H_pt)
  Ambient, _ = polynomial_ring(H_pt_loc, vcat(["e$i" for i in 1:n_vertices_G]))
  localized_cohomology, _ = quo(Ambient, ideal(Ambient, create_relations(gens(Ambient))))
  
  return GKMCohomology(
    H_pt,
    H_pt_loc,
    localized_cohomology,
    Dict{Edge,QQMPolyRingElem}(),
    QQMPolyRingElem[],
  )
end

