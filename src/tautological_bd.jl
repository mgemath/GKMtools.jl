export tautological_bd, univ_quotient_bd, wedge_product

@doc raw"""
    tautological_bd(::Type{GKM_graph}, k::Int, n::Int) -> GKM_vector_bundle

Return the tautological vector bundle of the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

# Example
```jldoctest
julia> S = tautological_bd(GKM_graph, 2, 4)
GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
12: (-1, 0, 0, 0), (0, -1, 0, 0)
13: (-1, 0, 0, 0), (0, 0, -1, 0)
14: (-1, 0, 0, 0), (0, 0, 0, -1)
23: (0, -1, 0, 0), (0, 0, -1, 0)
24: (0, -1, 0, 0), (0, 0, 0, -1)
34: (0, 0, -1, 0), (0, 0, 0, -1)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function tautological_bd(::Type{GKM_graph}, k::Int, n::Int)::GKM_vector_bundle
  return _tautological_bd(GKM_graph, [k, n-k], true)
end

@doc raw"""
    univ_quotient_bd(::Type{GKM_graph}, k::Int, n::Int) -> GKM_vector_bundle

Return the universal quotient bundle of the Grassmann variety of `k`-planes in the complex vector space of dimension `n`.

# Example
```jldoctest
julia> Q = univ_quotient_bd(GKM_graph, 2, 4)
GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
12: (0, 0, -1, 0), (0, 0, 0, -1)
13: (0, -1, 0, 0), (0, 0, 0, -1)
14: (0, -1, 0, 0), (0, 0, -1, 0)
23: (-1, 0, 0, 0), (0, 0, 0, -1)
24: (-1, 0, 0, 0), (0, 0, -1, 0)
34: (-1, 0, 0, 0), (0, -1, 0, 0)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function univ_quotient_bd(::Type{GKM_graph}, k::Int, n::Int)::GKM_vector_bundle
  return _tautological_bd(GKM_graph, [k, n-k], false)
end

function _tautological_bd(::Type{GKM_graph}, s::Vector{Int64}, tautological::Bool)::GKM_vector_bundle

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

  M_bd = M #free_module(ZZ, K[end])
  # matrix_morph = zero_matrix(ZZ, K[end], 2*K[end])
  GMtoM = ModuleHomomorphism(M, M_bd, [gens(M_bd)[i] for i in 1:K[end]])
  
  rank_of_bd = s[tautological ? 1 : 2] # s=[r, n-r], if tautologial is true, rank is r otherwise n-r
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(ZZ))}}(undef, nv, rank_of_bd)
  
  for v in 1:nv
    for r in 1:rank_of_bd
      weightMatrix[v, r] = -gens(M_bd)[d[v][r + s[1]*Int(!tautological)]]
    end
  end

  return vector_bundle(G, M_bd, GMtoM, weightMatrix; calculateConnection = true)
end


@doc raw"""
    wedge_product(V::GKM_vector_bundle, r::Int64) -> GKM_vector_bundle

Return the wedge product, or external product, `\wedge^r V`.

# Example
Let us compute the Plucker line bundle $l$ of the Grassmannian $G(2, 4)$.
```jldoctest
julia> S = tautological_bd(GKM_graph, 2, 4)
GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
12: (-1, 0, 0, 0), (0, -1, 0, 0)
13: (-1, 0, 0, 0), (0, 0, -1, 0)
14: (-1, 0, 0, 0), (0, 0, 0, -1)
23: (0, -1, 0, 0), (0, 0, -1, 0)
24: (0, -1, 0, 0), (0, 0, 0, -1)
34: (0, 0, -1, 0), (0, 0, 0, -1)

julia> l_minus_one = wedge_product(S, 2)
GKM vector bundle of rank 1 over GKM graph with 6 nodes and valency 4 with weights:
12: (-1, -1, 0, 0)
13: (-1, 0, -1, 0)
14: (-1, 0, 0, -1)
23: (0, -1, -1, 0)
24: (0, -1, 0, -1)
34: (0, 0, -1, -1)

julia> l = dual(l_minus_one)
GKM vector bundle of rank 1 over GKM graph with 6 nodes and valency 4 with weights:
12: (1, 1, 0, 0)
13: (1, 0, 1, 0)
14: (1, 0, 0, 1)
23: (0, 1, 1, 0)
24: (0, 1, 0, 1)
34: (0, 0, 1, 1)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function wedge_product(V::GKM_vector_bundle, r::Int64)::GKM_vector_bundle
  
  @req r<=rank(V) "r is greater than the rank"
  @req r>0 "r must be positive"

  G = V.gkm
  rank_w = binomial(rank(V), r)
  
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(_get_weight_type(G))}}(undef, n_vertices(G.g), rank_w)

  indices = collect(subsets(Set([i for i in 1:rank(V)]), r))

  for v in 1:n_vertices(G.g)
    for r in 1:rank_w
      weightMatrix[v, r] = sum(i -> V.w[v, i], indices[r])
    end
  end

  return vector_bundle(G, V.M, V.GMtoM, weightMatrix; calculateConnection = true)
end

@doc raw"""
    ^(V::GKM_vector_bundle, n::Number) -> GKM_vector_bundle

Return the tensor product `V^{\otimes n}`.

# Example
Let us compute the line bundle $l=\mathcal{O}(-4)$ of the Grassmannian $G(2, 4)$.
```jldoctest
julia> S = tautological_bd(GKM_graph, 2, 4)
GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
12: (-1, 0, 0, 0), (0, -1, 0, 0)
13: (-1, 0, 0, 0), (0, 0, -1, 0)
14: (-1, 0, 0, 0), (0, 0, 0, -1)
23: (0, -1, 0, 0), (0, 0, -1, 0)
24: (0, -1, 0, 0), (0, 0, 0, -1)
34: (0, 0, -1, 0), (0, 0, 0, -1)

julia> O_minus_one = wedge_product(S, 2)
GKM vector bundle of rank 1 over GKM graph with 6 nodes and valency 4 with weights:
12: (-1, -1, 0, 0)
13: (-1, 0, -1, 0)
14: (-1, 0, 0, -1)
23: (0, -1, -1, 0)
24: (0, -1, 0, -1)
34: (0, 0, -1, -1)

julia> l = O_minus_one^4
GKM vector bundle of rank 1 over GKM graph with 6 nodes and valency 4 with weights:
12: (-4, -4, 0, 0)
13: (-4, 0, -4, 0)
14: (-4, 0, 0, -4)
23: (0, -4, -4, 0)
24: (0, -4, 0, -4)
34: (0, 0, -4, -4)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function ^(V::GKM_vector_bundle, n::Number)::GKM_vector_bundle

  @req rank(V) == 1 "currently tensor product implemented only for line bundles"

  return vector_bundle(V.gkm, V.M, V.GMtoM, n*V.w; calculateConnection = true)
end