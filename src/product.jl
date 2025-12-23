@doc raw"""
    *(G1::AbstractGKM_graph, G2::AbstractGKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true) -> AbstractGKM_graph

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
function *(G1::AbstractGKM_graph, G2::AbstractGKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true)::AbstractGKM_graph

  if parent(_get_weight_type(G1)) == parent(_get_weight_type(G2))
    return _product(G1, G2; calculateCurveClasses, calculateConnection)
  end

  return _product(convert_weights(G1), convert_weights(G2); calculateCurveClasses, calculateConnection)
end

function _product(G1::AbstractGKM_graph, G2::AbstractGKM_graph; calculateCurveClasses::Bool=true, calculateConnection::Bool=true)::AbstractGKM_graph
  
  # @req G1.weightType == G2.weightType "GKM graphs must have the same weight type to be able to take their products"
  # @req base_ring(G1.M) == base_ring(G2.M) "GKM graphs must have the same character lattice base ring to be able to take their products"
  weightType = typeof(_get_weight_type(G1))
  baseRing = base_ring(G1.M)

  n1 = n_vertices(G1.g)
  n2 = n_vertices(G2.g)
  ne1 = n_edges(G1.g)
  ne2 = n_edges(G2.g)
  nv = n1 * n2

  g = Graph{Undirected}(nv)
  M = free_module(baseRing, rank(G1.M)+rank(G2.M)) # direct_sum(G1.M, G2.M)
  f1 = hom(G1.M, M, [gens(M)[i] for i in 1:rank(G1.M)])
  f2 = hom(G2.M, M, [gens(M)[i + rank(G1.M)] for i in 1:rank(G2.M)])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{weightType}}()
  labels = Vector{String}(undef, nv)

  # Store edge information for later connection building
  # Maps product edge to (source (:G1 or :G2), original_edge (on source), vertex_index (on the other Gi))
  edgeOrigin = Dict{Edge, Tuple{Symbol, Edge, Int64}}()  
  if calculateConnection
    # connection for product:
    con1 = get_connection(G1)
    con2 = get_connection(G2)
    if isnothing(con1) || isnothing(con2)
      calculateConnection = false
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

      # Store edge origin for connection building (both orientations)
      if calculateConnection
        edgeOrigin[E] = (:G1, e, _v)
        edgeOrigin[reverse(E)] = (:G1, reverse(e), _v) #TODO check if my update is okay here.
      end

      if calculateCurveClasses
        # here we assume that Oscar's direct sum preserves the order of generators
        edgeToGenIndex[E] = G1curveClasses.edgeToGenIndex[e] + (_v-1)*ne1
        edgeToGenIndex[reverse(E)] = edgeToGenIndex[E]
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

      # Store edge origin for connection building (both orientations)
      if calculateConnection
        edgeOrigin[E] = (:G2, e, v)
        edgeOrigin[reverse(E)] = (:G2, reverse(e), v)
      end

      if calculateCurveClasses
          # here we assume that Oscar's direct sum preserves the order of generators
          edgeToGenIndex[E] = offset + G2curveClasses.edgeToGenIndex[e] + (v-1)*ne2
          edgeToGenIndex[reverse(E)] = edgeToGenIndex[E]
      end
    end
  end

  # Build the GKM graph - need to handle standalone flags for non-compact case
  val1 = valency(G1)
  val2 = valency(G2)
  val_product = val1 + val2

  # Check if either graph is non-compact
  has_standalone = !is_compact(G1) || !is_compact(G2)

  # At least one graph has standalone flags - build manually
  weights_at_vertex = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{weightType}}}(undef, nv)
  flag_to_edge = Vector{Vector{Union{Nothing, Edge}}}(undef, nv)
  edge_to_flag_index = Dict{Edge, Int64}()
  for v in 1:n1
    for _v in 1:n2
      V = v + (_v-1)*n1
     # Flags at vertex V = (v, _v) consist of:
      # - First val1 flags from G1 (at the same _v coordinate)
      # - Next val2 flags from G2 (at the same v coordinate)
      weights_at_vertex[V] = Vector{AbstractAlgebra.Generic.FreeModuleElem{weightType}}(undef, val_product)
      flag_to_edge[V] = Vector{Union{Nothing, Edge}}(undef, val_product)
     # Copy flags from G1
      for i in 1:val1
        weights_at_vertex[V][i] = f1(G1.weights_at_vertex[v][i])
       # Check if this flag corresponds to an edge
        edge_G1 = G1.flag_to_edge[v][i]
        if !isnothing(edge_G1)
          # This is an edge flag - map to product
          @req src(edge_G1) == v "G1 has broken edge-flag relation. Check with isvalid."
          # v_other = src(edge_G1) == v ? dst(edge_G1) : src(edge_G1)
          v_other = dst(edge_G1)
          V_other = v_other + (_v-1)*n1
          E_product = Edge(V, V_other)
          flag_to_edge[V][i] = E_product
          edge_to_flag_index[E_product] = i
        else
          # Standalone flag
          flag_to_edge[V][i] = nothing
        end
      end
     # Copy flags from G2
      for i in 1:val2
        weights_at_vertex[V][val1 + i] = f2(G2.weights_at_vertex[_v][i])
       # Check if this flag corresponds to an edge
        edge_G2 = G2.flag_to_edge[_v][i]
        if !isnothing(edge_G2)
          # This is an edge flag - map to product
          @req src(edge_G2) == _v "G2 has broken edge-flag relation. Check with isvalid."
          # _v_other = src(edge_G2) == _v ? dst(edge_G2) : src(edge_G2)
          _v_other = dst(edge_G2)
          V_other = v + (_v_other-1)*n1
          E_product = Edge(V, V_other)
          flag_to_edge[V][val1 + i] = E_product
          edge_to_flag_index[E_product] = val1 + i
        else
          # Standalone flag
          flag_to_edge[V][val1 + i] = nothing
        end
      end
    end
  end
 # Build edge weight dict (already have this from earlier)
  for e in edges(g)
    W[reverse(e)] = -W[e]
  end
 # Create the GKM graph object
  GW_structure_consts = Dict{CurveClass_type, Array{Any, 3}}()
  res = AbstractGKM_graph(g, labels, M, weights_at_vertex, edge_to_flag_index, flag_to_edge, W,
                         nothing, nothing, nothing, GW_structure_consts, false)
  res.equivariantCohomology = _equivariant_cohomology_ring(res)

  if calculateCurveClasses
    dualConeRaySum, C, H2ToCN = _finish_GKM_H2(edgeLattice, H2, q, res, edgeToGenIndex)
    res.curveClasses = GKM_H2(res, edgeLattice, H2, edgeToGenIndex, q, dualConeRaySum, C, H2ToCN, nothing, nothing)
  end

  if calculateConnection
    # Build connection
    newCon = Dict{Edge, Vector{Int64}}()
    newA = Dict{Edge, Vector{ZZRingElem}}()

    # We need to build connections for both orientations of each edge
    all_oriented_edges = Vector{Edge}() # TODO: could iterate lazily instead.
    for e in edges(res.g)
      push!(all_oriented_edges, e)
      push!(all_oriented_edges, reverse(e))
    end

    # flags at vertex V=(v,_v) are [G1 flags at v, G2 flags at _v]
    for E in all_oriented_edges
      newCon[E] = Vector{Int64}(undef, val_product)
      newA[E] = Vector{ZZRingElem}(undef, val_product)
     (source, e_orig, param) = edgeOrigin[E]
      V1 = src(E)
      V2 = dst(E)
     # Get vertex coordinates in original graphs
      v1 = ((V1-1) % n1) + 1
      _v1 = div(V1-1, n1) + 1
      v2 = ((V2-1) % n1) + 1
      _v2 = div(V2-1, n1) + 1
     # For each flag at V1, determine where it connects to at V2
      for i_prod in 1:val_product
        if i_prod <= val1
          # This is a G1 flag (flag index i_orig in G1 at vertex v1)
          i_orig = i_prod
         if source == :G1
            # E is from G1: parallel connection
            # e_orig_oriented = src(e_orig) == v1 ? e_orig : reverse(e_orig)
            @req src(e_orig) == v1 "edge flag relation broken at G1"
            j_orig = con1.con[e_orig][i_orig]
            a_val = con1.a[e_orig][i_orig]
           # Flag j_orig in G1 at v2 maps to flag j_orig in product at V2=(v2,_v1)
            j_prod = j_orig
          else
            # E is from G2: orthogonal (trivial) connection
            # Flag i_prod at V1=(v1,_v1) stays as flag i_prod at V2=(v1,_v2)
            j_prod = i_prod
            a_val = ZZ(0)
          end
        else
          # This is a G2 flag (flag index i_orig in G2 at vertex _v1)
          i_orig = i_prod - val1
         if source == :G2
            # E is from G2: parallel connection
            # e_orig_oriented = src(e_orig) == _v1 ? e_orig : reverse(e_orig)
            @req src(e_orig) == _v1 "edge-flag relation is broken"
            j_orig = con2.con[e_orig][i_orig]
            a_val = con2.a[e_orig][i_orig]
           # Flag j_orig in G2 at _v2 maps to flag (val1 + j_orig) in product at V2=(v1,_v2)
            j_prod = val1 + j_orig
          else
            # E is from G1: orthogonal (trivial) connection
            # Flag i_prod at V1=(v1,_v1) stays as flag i_prod at V2=(v2,_v1)
            j_prod = i_prod
            a_val = ZZ(0)
          end
        end
        newCon[E][i_prod] = j_prod
        newA[E][i_prod] = a_val
      end
    end

    newConObj = GKM_connection(res, newCon, newA)
    set_connection!(res, newConObj)
  end

  return res
end
