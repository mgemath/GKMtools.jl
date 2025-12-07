export _tautological_bd

function _tautological_bd(::Type{GKM_graph}, s::Vector{Int64}, tautological::Bool)

  @req !isempty(s) "the vector of dimensions is empty"
  @req all(i -> s[i] > 0, eachindex(s)) "all dimensions must be positive"

  @req length(s) == 2 "this functionality is implemented for Grassmannians only"

  K::Vector{Int64} = [sum(s[1:i]) for i in 0:length(s)]
  d::Dict{Int64, NTuple{K[end], Int64}} = Dict{Int64, NTuple{K[end], Int64}}()
  index::Int64 = 1
  
  for c in Combinatorics.permutations(1:K[end])
    if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:length(s))
      d[index] = (c...,)
      index += 1
    end
  end
      
  nv::Int64 = length(keys(d))
  g = Graph{Undirected}(nv)

  M = free_module(ZZ, K[end])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()

  get_root::Dict{Edge, Vector{Int64}} = Dict{Edge, Vector{Int64}}()

  for v in 1:nv
    for w in (v+1):nv
      dif = _T_link(d[v], d[w], K)
        
      isempty(dif) && continue

      add_edge!(g, v, w)
      W[Edge(w, v)] = gens(M)[dif[2]] - gens(M)[dif[1]]

      get_root[Edge(w, v)] = dif
    end
  end

  ## construct connection
  a::Dict{Tuple{Edge, Edge}, ZZRingElem} = Dict{Tuple{Edge, Edge}, ZZRingElem}()

  for _v in vertices(g)
    for _w in all_neighbors(g, _v)

      alpha = get_root[_w < _v ? Edge(_v, _w) : Edge(_w, _v)]


      for _u in all_neighbors(g, _v)
        # e = (_v, _w)
        # e'= (_v, _u)

        beta = get_root[_u < _v ? Edge(_v, _u) : Edge(_u, _v)]

        if alpha == beta
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ(2)
        elseif alpha[1] == beta[1] || alpha[2] == beta[2]
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ((_v < _w ? -1 : 1)*(_v < _u ? -1 : 1))
        elseif alpha[1] == beta[2] || alpha[2] == beta[1]
          a[(Edge(_v, _w), Edge(_v, _u))] = -ZZ((_v < _w ? -1 : 1)*(_v < _u ? -1 : 1))
        else
          a[(Edge(_v, _w), Edge(_v, _u))] = ZZ(0)
        end
        
      end
    end
  end

  labels = [prod(i -> d[n][i]>9 ? "|$(d[n][i])|" : "$(d[n][i])", 1:K[end-1]) for n in 1:nv] 

  G = gkm_graph(g, labels, M, W)


  con = build_GKM_connection(G, a)
  set_connection!(G, con)

  M_bd = free_module(ZZ, 2*K[end])
  # matrix_morph = zero_matrix(ZZ, K[end], 2*K[end])
  GMtoM = ModuleHomomorphism(M, M_bd, [gens(M_bd)[i] for i in 1:K[end]])
  
  rank_of_bd = s[tautological ? 1 : 2] # s=[r, n-r], if tautologial is true, rank is r otherwise n-r
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(ZZ))}}(undef, nv, rank_of_bd)
  
  for v in 1:nv
    for r in 1:rank_of_bd
      ## add the weights of the basis
      weightMatrix[v, r] = sum(l -> gens(M_bd)[d[v][l]], 1:s[1])

      ## add the weights of the fiber
      weightMatrix[v, r] = weightMatrix[v, r] + gens(M_bd)[d[v][r + s[1]*Int(!tautological)] + rank_torus(G)]
    end
  end

  return vector_bundle(G, M_bd, GMtoM, weightMatrix; calculateConnection = true)
end