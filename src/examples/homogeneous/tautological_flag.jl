@doc raw"""
    tautological_and_univ_bd(::Type{GKMGraph}, k::Integer, n::Integer)

Return the tautological bundle and universal quotient bundle on the
Grassmannian `Gr(k,n)`. The base graph is built once, with its homogeneous
connection computed in place, and the stored flag representative at every
vertex determines the two fibre weight spaces.

# Examples
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKMGraph, 2, 4);

julia> S
GKM vector bundle of rank 2 over 6 vertices with weights:
12: (-1, 0, 0, 0), (0, -1, 0, 0)
13: (-1, 0, 0, 0), (0, 0, -1, 0)
14: (-1, 0, 0, 0), (0, 0, 0, -1)
23: (0, -1, 0, 0), (0, 0, -1, 0)
24: (0, -1, 0, 0), (0, 0, 0, -1)
34: (0, 0, -1, 0), (0, 0, 0, -1)

julia> Q
GKM vector bundle of rank 2 over 6 vertices with weights:
12: (0, 0, -1, 0), (0, 0, 0, -1)
13: (0, -1, 0, 0), (0, 0, 0, -1)
14: (0, -1, 0, 0), (0, 0, -1, 0)
23: (-1, 0, 0, 0), (0, 0, 0, -1)
24: (-1, 0, 0, 0), (0, 0, -1, 0)
34: (-1, 0, 0, 0), (0, -1, 0, 0)
```
"""
function tautological_and_univ_bd(
  ::Type{GKMGraph}, k::Integer, n::Integer,
)
  @req 0 < k < n "require 0 < k < n"

  G = grassmannian(GKMGraph, k, n)
  M = lattice(G)
  basis = gens(M)
  GMtoM = hom(M, M, basis)
  tautological_weights = Matrix{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(
    undef, num_vertices(G), k,
  )
  quotient_weights = Matrix{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(
    undef, num_vertices(G), n - k,
  )

  for v in vertices(G)
    representative = vertices_structure(G)[v].flag
    for i in 1:k
      tautological_weights[v, i] = -basis[representative[i]]
    end
    for i in 1:(n - k)
      quotient_weights[v, i] = -basis[representative[k + i]]
    end
  end

  return (
    vector_bundle(G, M, GMtoM, tautological_weights),
    vector_bundle(G, M, GMtoM, quotient_weights),
  )
end
