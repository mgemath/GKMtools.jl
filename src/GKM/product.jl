@doc raw"""
    *(G1::GKM_graph, G2::GKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true) -> GKM_graph

It constructs the product of two GKM graphs.

# Examples
```jldoctest
julia> G = generalized_gkm_flag(root_system(:A, 1))
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

julia> G*G
GKM graph with 4 nodes, valency 2 and axial function:
s1,id -> id,id => (-1, 1, 0, 0)
id,s1 -> id,id => (0, 0, -1, 1)
s1,s1 -> s1,id => (0, 0, -1, 1)
s1,s1 -> id,s1 => (-1, 1, 0, 0)

```

!!! note
    The character group is of type free ``\mathbb{Q}``-module if this holds for one of the two GKM graphs.
"""
function *(G1::GKM_graph, G2::GKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true)::GKM_graph

  if _get_weight_type(G1) == _get_weight_type(G2)
    return _product(G1, G2; calculateCurveClasses, calculateConnection)
  end

  return _product(convert_weights(G1), convert_weights(G2); calculateCurveClasses, calculateConnection)
end

function _product(G1::GKM_graph, G2::GKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true)::GKM_graph
  
  # @req G1.weightType == G2.weightType "GKM graphs must have the same weight type to be able to take their products"
  # @req base_ring(G1.M) == base_ring(G2.M) "GKM graphs must have the same character lattice base ring to be able to take their products"
  weightType = typeof(_get_weight_type(G1))
  baseRing = base_ring(G1.M)

  n1 = n_vertices(G1.g)
  n2 = n_vertices(G2.g)
  ne1 = n_edges(G1.g)
  ne2 = n_edges(G2.g)
  nv = n1 * n_vertices(G2.g)

  g = Graph{Undirected}(nv)
  M = free_module(baseRing, rank(G1.M)+rank(G2.M)) # direct_sum(G1.M, G2.M)
  f1 = hom(G1.M, M, [gens(M)[i] for i in 1:rank(G1.M)])
  f2 = hom(G2.M, M, [gens(M)[i + rank(G1.M)] for i in 1:rank(G2.M)])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{weightType}}()
  labels = Vector{String}(undef, nv)

  if calculateConnection
    # connection for product:
    con1 = get_connection(G1)
    con2 = get_connection(G2)
    if isnothing(con1) || isnothing(con2)
      calculateConnection = false
    else
      newCon = Dict{Tuple{Edge, Edge}, Edge}()
    end
  end

  # only needed if calculateCurveClasses==true
  if calculateCurveClasses

    G1curveClasses = GKM_second_homology(G1)
    G2curveClasses = GKM_second_homology(G2)

    edgeLattice1 = G1curveClasses.edgeLattice
    edgeLattice2 = G2curveClasses.edgeLattice
    G1H2 = G1curveClasses.H2
    G2H2 = G2curveClasses.H2
    q1 = G1curveClasses.quotientMap
    q2 = G2curveClasses.quotientMap

    # first all G1-edges for each vertex of G2, then vice versa.
    edgeLattice, _, _ = direct_sum(vcat(repeat([edgeLattice1], n2), repeat([edgeLattice2], n1)))
    H2, _, _ = direct_sum(G1H2, G2H2)
    qMatrix = vcat(repeat([q1 Int(0)], n2), repeat([Int(0) q2], n1))
    q = ModuleHomomorphism(edgeLattice, H2, qMatrix) # direct sum of morphisms
    
    edgeToGenIndex = Dict{Edge, Int64}()
  end

  # labels
  for v in 1:n1
    for _v in 1:n2
      V = v + (_v-1)*n1
      labels[V] = G1.labels[v]*","*G2.labels[_v]
    end
  end

  # edges from G1
  for e in edges(G1.g)
    v, w = src(e), dst(e)
    for _v in 1:n2
      V1 = v + (_v-1)*n1
      V2 = w + (_v-1)*n1
      add_edge!(g, V1, V2)
      E = Edge(V1, V2)
      W[E] = f1(G1.w[e])        
      if calculateCurveClasses
        # here we assume that Oscar's direct sum preserves the order of generators
        edgeToGenIndex[E] = G1curveClasses.edgeToGenIndex[e] + (_v-1)*ne1
        edgeToGenIndex[reverse(E)] = edgeToGenIndex[E]
      end
      if calculateConnection
        #first, copy old connection:
        for u in all_neighbors(G1.g, v)
          ei = Edge(v, u)
          epi = con1.con[(e, ei)]
          U1 = u + (_v-1)*n1
          U2 = dst(epi) + (_v-1)*n1
          Ei = Edge(V1, U1)
          Epi = Edge(V2, U2)
          newCon[(E, Ei)] = Epi
          newCon[(reverse(E), Epi)] = Ei
        end
        #second, add trivial connection in normal direction:
        for _u in all_neighbors(G2.g, _v)
          U1 = v + (_u-1)*n1
          U2 = w + (_u-1)*n1
          Ei = Edge(V1, U1)
          Epi = Edge(V2, U2)
          newCon[(E, Ei)] = Epi
          newCon[(reverse(E), Epi)] = Ei
        end
      end
    end
  end

  # edges from G2
  offset = ne1 * n2
  for e in edges(G2.g)
    _v, _w = src(e), dst(e)
    for v in 1:n1
      V1 = v + (_v-1)*n1
      V2 = v + (_w-1)*n1
      add_edge!(g, V1, V2)
      E = Edge(V1, V2)
      W[E] = f2(G2.w[e])
      if calculateCurveClasses
          # here we assume that Oscar's direct sum preserves the order of generators
          edgeToGenIndex[E] = offset + G2curveClasses.edgeToGenIndex[e] + (v-1)*ne2
          edgeToGenIndex[reverse(E)] = edgeToGenIndex[E]
      end
      if calculateConnection
        #first, copy old connection:
        for _u in all_neighbors(G2.g, _v)
          ei = Edge(_v, _u)
          epi = con2.con[(e, ei)]
          U1 = v + (_u-1)*n1
          U2 = v + (dst(epi)-1)*n1
          Ei = Edge(V1, U1)
          Epi = Edge(V2, U2)
          newCon[(E, Ei)] = Epi
          newCon[(reverse(E), Epi)] = Ei
        end
        #second, add trivial connection in normal direction:
        for u in all_neighbors(G1.g, v)
          U1 = u + (_v-1)*n1
          U2 = u + (_w-1)*n1
          Ei = Edge(V1, U1)
          Epi = Edge(V2, U2)
          newCon[(E, Ei)] = Epi
          newCon[(reverse(E), Epi)] = Ei
        end
      end
    end
  end

  res = gkm_graph(g, labels, M, W)

  if calculateCurveClasses
    dualConeRaySum, C, H2ToCN = _finish_GKM_H2(edgeLattice, H2, q, res, edgeToGenIndex)
    res.H2 = GKM_H2(edgeLattice, H2, edgeToGenIndex, q, dualConeRaySum, C, H2ToCN, nothing, nothing)
  end
  if calculateConnection
    newConObj = build_gkm_connection(res, newCon)
    set_connection!(res, newConObj)
  end

  return res
end
