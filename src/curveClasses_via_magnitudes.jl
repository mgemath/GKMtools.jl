@doc raw"""
    _primitive_edge_integral_map(
      G::AbstractGKM_graph,
      class::FreeModElem{QQMPolyRingElem}=first_chern_class(G)
    )

Return a surjective homomorphism from `free_module(ZZ, n_edges(G.g))` to
`free_module(ZZ, 1)` obtained by integrating `class` over the GKM graph edges.
The `i`-th source generator corresponds to the `i`-th edge in
`collect(edges(G.g))`.

The edge integrals are divided by their greatest common divisor. An error is
thrown when there are no edges or when all edge integrals vanish, since in
those cases no such surjective homomorphism exists.
"""
function _primitive_edge_integral_map(
  G::AbstractGKM_graph,
  class::FreeModElem{QQMPolyRingElem}=first_chern_class(G)
)::AbstractAlgebra.Generic.ModuleHomomorphism{ZZRingElem}
  cohomologyRing = G.equivariantCohomology
  @req parent(class) === cohomologyRing.cohomRing "The cohomology class must belong to the equivariant cohomology of G"

  edgeList = collect(edges(G.g))
  @req !isempty(edgeList) "The GKM graph must have at least one edge"

  edgeIntegrals = [_integrate_edge_to_int(class, G, e) for e in edgeList]
  commonDivisor = foldl(gcd, edgeIntegrals; init=ZZ(0))
  # If this is called for computing H2, the reason for the error above is that `class` integrates to zero over each edge.
  # For GKM graphs of compact Hamiltonian GKM spaces this should never happen when `class` is the first Chern class.
  @req !iszero(commonDivisor) "All edge integrals vanish, so the resulting homomorphism cannot be surjective"

  primitiveIntegrals = [divexact(value, commonDivisor) for value in edgeIntegrals]
  source = free_module(ZZ, length(edgeList))
  target = free_module(ZZ, 1)
  targetGen = gens(target)[1]
  images = [value * targetGen for value in primitiveIntegrals]

  return ModuleHomomorphism(source, target, images)
end
