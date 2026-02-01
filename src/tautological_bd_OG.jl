export tautological_and_univ_bd_og

# @doc raw"""
#     tautological_and_univ_bd(::Type{GKM_graph}, k::Int, n::Int) -> GKM_vector_bundle{ZZRingElem}

# Return a pair `(S, Q)` where `S` is the tautological vector bundle, and `Q` is the universal quotient bundle of the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

# # Example
# ```jldoctest
# julia> S, Q = tautological_and_univ_bd(GKM_graph, 2, 4);

# julia> S
# GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
# 12: (-1, 0, 0, 0), (0, -1, 0, 0)
# 13: (-1, 0, 0, 0), (0, 0, -1, 0)
# 14: (-1, 0, 0, 0), (0, 0, 0, -1)
# 23: (0, -1, 0, 0), (0, 0, -1, 0)
# 24: (0, -1, 0, 0), (0, 0, 0, -1)
# 34: (0, 0, -1, 0), (0, 0, 0, -1)

# julia> Q
# GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
# 12: (0, 0, -1, 0), (0, 0, 0, -1)
# 13: (0, -1, 0, 0), (0, 0, 0, -1)
# 14: (0, -1, 0, 0), (0, 0, -1, 0)
# 23: (-1, 0, 0, 0), (0, 0, 0, -1)
# 24: (-1, 0, 0, 0), (0, 0, -1, 0)
# 34: (-1, 0, 0, 0), (0, -1, 0, 0)
# ```
# !!! warning
#     All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

# """
function tautological_and_univ_bd_og(::Type{GKM_graph}, k::Int, n::Int)::Tuple{GKM_vector_bundle{ZZRingElem}, GKM_vector_bundle{ZZRingElem}}
  @req k > 0 "k must be positive"
  @req n > k "n must be greater than k"
  @req iseven(n) "n must be even for OG Grassmannians"
  return _tautological_and_univ_bd_og(GKM_graph, [k, n-k])
end
function _tautological_and_univ_bd_og(::Type{GKM_graph}, s::Vector{Int64})::Tuple{GKM_vector_bundle{ZZRingElem}, GKM_vector_bundle{ZZRingElem}}

  @req !isempty(s) "the vector of dimensions is empty"
  @req all(i -> s[i] > 0, eachindex(s)) "all dimensions must be positive"

  @req length(s) == 2 "this functionality is implemented for Grassmannians only"

  K::Vector{Int64} = [sum(s[1:i]) for i in 0:length(s)]
  d::Dict{Int64, NTuple{K[end], Int64}} = Dict{Int64, NTuple{K[end], Int64}}()
  index::Int64 = 1
  
  for c in Combinatorics.permutations(1:K[end])
    if all(i-> issorted(c[(K[i]+1):K[i+1]]), 1:length(s))
      # check isotropic condition
      _is_isotropic(c, s) || continue

      # _is_symplectic(c, s) || continue

      d[index] = (c...,)
      index += 1
    end
  end
  println("Number of fixed points: $(length(keys(d)))")
  println(d)
  nv::Int64 = length(keys(d))
  g = Graph{Undirected}(nv)

  M = free_module(ZZ, K[end])
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()

  get_root::Dict{Edge, Vector{Int64}} = Dict{Edge, Vector{Int64}}()

  for v in 1:nv
    for w in (v+1):nv

      ###############
      # compute the difference
      # you cannot have a curve between flag <e1> and flag <e8> in C8, but you can have it between <e1> and <e7>, because e8 is the dual of e1
      # you can have a curve between <e1, e2> and <e5, e6>, because none of them are dual to each other
      node_v = d[v]
      node_w = d[w]
      found = 0
      for i in 1:s[1]
        if any(j -> node_v[i] + node_w[j] == K[end] + 1, 1:s[1])
          found += 1
        end
      end
      isodd(found) && continue  
      ###############
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

  # con = build_GKM_connection(G, a)
  # set_connection!(G, con)

  M_bd = M #free_module(ZZ, K[end])
  # matrix_morph = zero_matrix(ZZ, K[end], 2*K[end])
  GMtoM = ModuleHomomorphism(M, M_bd, [gens(M_bd)[i] for i in 1:K[end]])
  
  # rank_of_bd = s[tautological ? 1 : 2] # s=[r, n-r], if tautologial is true, rank is r otherwise n-r
  weightMatrix_S = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(ZZ))}}(undef, nv, s[1])
  weightMatrix_Q = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(ZZ))}}(undef, nv, s[2])
  
  for v in 1:nv
    for r in 1:s[1]
      weightMatrix_S[v, r] = -gens(M_bd)[d[v][r]]
    end
    for r in 1:s[2]
      weightMatrix_Q[v, r] = -gens(M_bd)[d[v][r + s[1]]]
    end
    # for r in 1:rank_of_bd
    #   weightMatrix[v, r] = -gens(M_bd)[d[v][r + s[1]*Int(!tautological)]]
    # end
  end

  S = vector_bundle(G, M_bd, GMtoM, weightMatrix_S; calculateConnection = true)
  Q = vector_bundle(G, M_bd, GMtoM, weightMatrix_Q; calculateConnection = true)

  return (S, Q)
end

function _is_isotropic(c, s)::Bool
  n = sum(s)
  for i in 1:s[1]
    any(j -> c[i] + c[j] == n+1, 1:s[1]) && return false
  end
  return true
  
end

function _is_symplectic(c, s)::Bool
  n = sum(s)
  ans = true
  for i in 1:s[1]
    any(j -> c[i] + c[j] == n+1, 1:s[1]) && continue
    ans = false
    break
  end
  return ans
end