function _evaluate_weight_functional(
  coefficients::Vector{R},
  weight::AbstractAlgebra.Generic.FreeModuleElem{R}
)::R where R <: GKM_weight_type
  value = zero(coefficients[1])
  for i in eachindex(coefficients)
    value += coefficients[i] * weight[i]
  end
  return value
end

@doc raw"""
    _two_independent_projection(G::AbstractGKM_graph)

Return a deterministic surjective module homomorphism from `G.M` to a free
module of rank two over the same coefficient ring such that the images of
every two flag weights at a vertex remain linearly independent.

The input graph must be 2-independent and its weight module must have rank at
least two. For integer weights, the map is surjective as a homomorphism of
integer lattices, not just after tensoring with the rationals.
"""
function _two_independent_projection(
  G::AbstractGKM_graph{R}
)::AbstractAlgebra.Generic.ModuleHomomorphism{R} where R <: GKM_weight_type
  r = rank_torus(G)
  @req r >= 2 "The weight module must have rank at least two"
  @req all(w -> !iszero(w), Iterators.flatten(G.weights_at_vertex)) "Flag weights must be nonzero"
  @req valency(G) < 2 || is2_indep(G) "The GKM graph must be 2-independent"

  coefficientRing = base_ring(G.M)

  # Find a functional a which is nonzero on every flag weight. For each
  # nonzero weight w, a(t)(w), with a(t) = (1, t, ..., t^(r-1)), is a
  # nonzero polynomial. Consequently only finitely many integer t are bad.
  firstCoefficients = Vector{R}()
  firstValues = Vector{Vector{R}}()
  t = 0
  while true
    firstCoefficients = [coefficientRing(t)^(i - 1) for i in 1:r]
    firstValues = [
      [_evaluate_weight_functional(firstCoefficients, w) for w in weights]
      for weights in G.weights_at_vertex
    ]
    all(value -> !iszero(value), Iterators.flatten(firstValues)) && break
    t += 1
  end

  # Once a is fixed, a second functional b must avoid one linear hyperplane
  # for every pair (w_i, w_j):
  #
  #   a(w_i)b(w_j) - a(w_j)b(w_i) != 0.
  #
  # Searching b(s) = (0, 1, s, ..., s^(r-2)) terminates deterministically.
  # Indeed, an identically vanishing polynomial would force the vector
  # a(w_i)w_j - a(w_j)w_i to be supported in its first coordinate. It is
  # also annihilated by a, whose first coordinate is one, so it would be
  # zero, contradicting the independence of w_i and w_j.
  secondCoefficients = Vector{R}()
  s = 0
  while true
    secondCoefficients = [
      i == 1 ? zero(coefficientRing) : coefficientRing(s)^(i - 2)
      for i in 1:r
    ]
    secondValues = [
      [_evaluate_weight_functional(secondCoefficients, w) for w in weights]
      for weights in G.weights_at_vertex
    ]

    preservesIndependence = true
    for v in eachindex(G.weights_at_vertex)
      weights = G.weights_at_vertex[v]
      for i in 1:length(weights)
        for j in (i + 1):length(weights)
          determinant = firstValues[v][i] * secondValues[v][j]
          determinant -= firstValues[v][j] * secondValues[v][i]
          if iszero(determinant)
            preservesIndependence = false
            break
          end
        end
        preservesIndependence || break
      end
      preservesIndependence || break
    end

    preservesIndependence && break
    s += 1
  end

  target = free_module(coefficientRing, 2)
  targetGens = gens(target)
  images = [
    firstCoefficients[i] * targetGens[1] + secondCoefficients[i] * targetGens[2]
    for i in 1:r
  ]
  return ModuleHomomorphism(G.M, target, images)
end
