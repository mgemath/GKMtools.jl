@doc raw"""
    Seidel_space(G::AbstractGKM_graph, w::AbstractAlgebra.Generic.FreeModuleElem{R}; basePoint::Int64 = 1) -> AbstractGKM_graph

Construct the Seidel space associated to the GKM graph `G` (representing the GKM variety $X$) and the map $\iota:\mathbb{C}^\times\rightarrow T$ given by the element $w\in\mathfrak{t}$.

# Optional arguments:
 - `basePoint::Int64`: This is the vertex of `G` so that in the internal presentation of the curve classes of $S_X$, the curve class of the section of $S_X\rightarrow \mathbb{P}^1$
    associated to the vertex `basePoint` is represented by `(0,...,0,1)`.
    The first entries correspond to curve classes of $X$. The last is the degree of the curve class projected to $\mathbb{P}^1$.

# Examples
```jldoctest Seidel_space
julia> G = projective_space(GKM_graph, 2);

julia> S = Seidel_space(G, gens(G.M)[1])
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, 1)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 0)

julia> S = Seidel_space(G, gens(G.M)[1] + 7*gens(G.M)[2])
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, -6)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 7)

julia> print_curve_classes(S)
[2]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [2]_0: (1, 0), Chern number: 3
[1]_inf -> [1]_0: (0, 1), Chern number: -3
[2]_inf -> [2]_0: (6, 1), Chern number: 15
[2]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [3]_0: (-1, 1), Chern number: -6
[3]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [2]_inf: (1, 0), Chern number: 3

```
Using a different base point does not change the resulting GKM graph but gives a different internal presentation of the curve classes.
```jldoctest Seidel_space
julia> S = Seidel_space(G, gens(G.M)[1] + 7*gens(G.M)[2]; basePoint=2)
GKM graph with 6 nodes, valency 3 and axial function:
[2]_0 -> [1]_0 => (-1, 1, 0, 0)
[3]_0 -> [1]_0 => (-1, 0, 1, 0)
[3]_0 -> [2]_0 => (0, -1, 1, 0)
[1]_inf -> [1]_0 => (0, 0, 0, -1)
[2]_inf -> [2]_0 => (0, 0, 0, -1)
[2]_inf -> [1]_inf => (-1, 1, 0, -6)
[3]_inf -> [3]_0 => (0, 0, 0, -1)
[3]_inf -> [1]_inf => (-1, 0, 1, 1)
[3]_inf -> [2]_inf => (0, -1, 1, 7)

julia> print_curve_classes(S)
[2]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [1]_0: (1, 0), Chern number: 3
[3]_0 -> [2]_0: (1, 0), Chern number: 3
[1]_inf -> [1]_0: (-6, 1), Chern number: -3
[2]_inf -> [2]_0: (0, 1), Chern number: 15
[2]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [3]_0: (-7, 1), Chern number: -6
[3]_inf -> [1]_inf: (1, 0), Chern number: 3
[3]_inf -> [2]_inf: (1, 0), Chern number: 3
```
"""
function Seidel_space(
  G::AbstractGKM_graph,
  weight::AbstractAlgebra.Generic.FreeModuleElem{R};
  basePoint::Int64 = 1
) where R<:GKM_weight_type

  @req parent(weight) === G.M "The weight does not belong to the right character lattice"
  @req is_connected(G.g) "The GKM graph should be connected."

  nv = n_vertices(G.g)
  ne = n_edges(G.g)
  r = rank_torus(G)
  val = valency(G)
  labels = G.labels
  
  # create labels for Seidel space (vertices 1...n over 0, n+1...2n over inf)
  Slabels = Vector{String}()
  sizehint!(Slabels, 2*nv)
  for l in labels
    push!(Slabels, "[" * l * "]_0")
  end
  for l in labels
    push!(Slabels, "[" * l * "]_inf")
  end

  # weights:
  SM = free_module(base_ring(G.M), r+1)
  z = gens(SM)[r+1]
  iM = ModuleHomomorphism(G.M, SM, [gens(SM)[i] for i in 1:r])
  SW = Vector{Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}}(undef, 2*nv)

  # copy weights of flags at all vertices
  for v in 1:nv
    w_at_v = G.weights_at_vertex[v]
    SW[v] = [iM(w_at_v[i]) for i in 1:val]   
    SW[v+nv] = [iM(w_at_v[i]) - sum([w_at_v[i][j] * weight[j] for j in 1:r]) * z for i in 1:val]
  end
  
  # new GKM graph:
  SG = flags_only_gkm_graph(Slabels, SM, SW; checkLabels=false)

  # add edges over zero and infinity.
  for e in edges(G.g)
    s = src(e)
    d = dst(e)
    ind_at_s = G.edge_to_flag_index[e]
    ind_at_d = G.edge_to_flag_index[reverse(e)]
    connect_flags!(SG, s, d, ind_at_s, ind_at_d)
    connect_flags!(SG, s+nv, d+nv, ind_at_s, ind_at_d)
  end

  # add horizontal edges
  for v in 1:nv
    add_edge!(SG, v, v+nv, z)
  end

  # infer GKM connection to Seidel space if possible
  con = get_connection(G)
  if !isnothing(con)
    SconDict = Dict{Edge, Vector{Int64}}()

    # Copy original connection for each original edge
    for e_unoriented in edges(G.g)
      for e in [e_unoriented, reverse(e_unoriented)]
        s = src(e)
        d = dst(e)
        e_0 = e
        e_inf = Edge(s + nv, d + nv)
        # Copy original connection for first val edges, and parallel connection for horizontal edge.
        SconDict[e_0] = [i <= val ? con.con[e][i] : val + 1 for i in 1:val+1]
        SconDict[e_inf] = [i <= val ? con.con[e][i] : val + 1 for i in 1:val+1]
      end
    end

    # Create connection over horizontal edges
    for v in 1:nv
      e = Edge(v, nv+v)
      SconDict[e] = [i for i in 1:val+1]
      SconDict[reverse(e)] = [i for i in 1:val+1]
    end

    Scon = build_GKM_connection(SG, SconDict)
    set_connection!(SG, Scon)
  end

  # build curve classes:

  GCC = GKM_second_homology(G)
  Gq = GCC.quotientMap
  rkGH2 = length(gens(GCC.H2))
  horizontalEdgeLattice = free_module(ZZ, nv)
  # first edges over 0, then over inf, then horizontal ones.
  SedgeLattice, _, _ = direct_sum(vcat(repeat([GCC.edgeLattice], 2), [horizontalEdgeLattice]))
  extraFactor = free_module(ZZ, 1)
  SH2, _, projs = direct_sum(GCC.H2, extraFactor)

  #build map from horizontalEdgeLattice to GCC.edgeLattice and extraFactor
  horToExtra = ModuleHomomorphism(horizontalEdgeLattice, extraFactor, [gens(extraFactor)[1] for i in 1:nv])
  # calculate H2 shifts of horizontal edges
  shiftDict = Dict{Int64, CurveClass_type}()
  shiftDict[1] = zero(GCC.H2)
  #first, assume basePoint = 1
  while length(keys(shiftDict)) < nv # this is why we checked connectedness of G first.
    for v1 in collect(keys(shiftDict))
      for v2 in all_neighbors(G.g, v1)
        haskey(shiftDict, v2) && continue
        e = Edge(v1,v2)
        we = _w(G, e)
        wz = sum([we[i] * weight[i] for i in 1:r])
        shiftDict[v2] = shiftDict[v1] -wz * curve_class(G, e)
      end
    end
  end
  # Now shift so that the correct base point has zero vertical class.
  horToGH2 = ModuleHomomorphism(horizontalEdgeLattice, GCC.H2, [shiftDict[i] - shiftDict[basePoint] for i in 1:nv])
  qMatrix = [Gq Int(0); Gq Int(0); horToGH2 horToExtra ]
  Sq = ModuleHomomorphism(SedgeLattice, SH2, qMatrix)
  SedgeToGenIndex = Dict{Edge, Int64}()
  for e in edges(G.g)
    indE = GCC.edgeToGenIndex[e]
    SedgeToGenIndex[e] = indE
    SedgeToGenIndex[reverse(e)] = indE
    SedgeToGenIndex[Edge(src(e)+nv, dst(e)+nv)] = indE + ne
    SedgeToGenIndex[Edge(dst(e)+nv, src(e)+nv)] = indE + ne
  end
  for v in 1:nv
    SedgeToGenIndex[Edge(v, v+nv)] = 2*ne + v
    SedgeToGenIndex[Edge(v+nv, v)] = 2*ne + v
  end
  dualConeRaySum, C, H2ToCN = _finish_GKM_H2(SedgeLattice, SH2, Sq, SG, SedgeToGenIndex)
  SG.curveClasses = GKM_H2(SG, SedgeLattice, SH2, SedgeToGenIndex, Sq, dualConeRaySum, C, H2ToCN, nothing, projs[1])

  return SG
end

function _SeidelSectionCount(SG::AbstractGKM_graph)
  @req divides(n_vertices(SG.g), 2)[1] "SG is not a Seidel space!"
  H2 = GKM_second_homology(SG)
  rkH2 = length(gens(H2.H2)) - 1

  if !isnothing(H2.sectionCount)
    return H2.sectionCount
  end

  # generate sectional multiplicity as ZZ-module homorphism from H2 to ZZ
  ZZasModule = free_module(ZZ, 1)
  H2.sectionCount = ModuleHomomorphism(H2.H2, ZZasModule, vcat([zero(ZZasModule) for i in 1:rkH2], [gens(ZZasModule)[1]]))
  return H2.sectionCount
end

# Warning: Don't use this on anything that is not the output of Seidel_space(...).
function _effectiveSectionClassesWithChernNumber(
  SG::AbstractGKM_graph,
  chernNumber::ZZRingElem;
)

  H2 = GKM_second_homology(SG)
  secCt = _SeidelSectionCount(SG)
  cN = H2.chernNumber

  ZZ2, _, _ = direct_sum([codomain(cN), codomain(secCt)])
  q = ModuleHomomorphism(H2.H2, ZZ2, hcat(matrix(cN), matrix(secCt)))

  success, e0 = has_preimage_with_preimage(q, chernNumber * gens(ZZ2)[1] + gens(ZZ2)[2])
  !success && return Vector{}()

  # e0 = nothing
  # try
  #   e0 = preimage(q, chernNumber * gens(ZZ2)[1] + gens(ZZ2)[2])
  # catch err
  #   if isa(err, ArgumentError)
  #     # In this case, no integer combination of edge curve classes has this chern number.
  #     return Vector{}()
  #   else
  #     rethrow(err)
  #   end
  # end

  K, k = kernel(q)
  rk = rank(K)
  mk = transpose(matrix(k)) # columns are images of generators

  dualRays = rays(H2.dualCone)
  nDualRays = length(dualRays)
  rkH2 = length(gens(H2.H2))
  rayMatrix = QQMatrix(nDualRays, rkH2)
  for i in 1:nDualRays
    for j in 1:rkH2
      rayMatrix[i, j] = dualRays[i][j]
    end
  end
  
  Re0 = rayMatrix*[e0[i] for i in 1:rkH2]
  P = polyhedron(-rayMatrix*mk, Re0)

  ptsIterator = (e0 + k(K([v[i] for i in 1:rk])) for v in lattice_points(P))
  return ptsIterator
end
