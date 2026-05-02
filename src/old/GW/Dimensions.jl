export dim_moduli_space

function dim_moduli_space(G::AbstractGKM_graph, beta::CurveClass_type, g::Int64, m::Int64)
  @req g >= 0 "Genus must be non-negative"
  @req m >= 0 "Number of marked points must be non-negative"

  return (1-g) * (valency(G) - 3) + m + chern_number(G, beta)
end