@doc raw"""
    line_bundle(G::AbstractGKM_graph, M::AbstractAlgebra.Generic.FreeModule{R}, GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}, weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}) -> GKM_vector_bundle

Return the equivariant line bundle over the GKM graph `G` whose weights over the fixed points are given by `weights`.

# Arguments
- `g::G::AbstractGKM_graph`: A GKM graph
- `M::AbstractAlgebra.Generic.FreeModule{R}`: The weight lattice of the torus acting on the line bundle.
    This is often bigger than the torus acting on `G`, for example when there is an extra scaling-action on the fibres.
- `GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}`: The inclusion of `G.M` (the weight lattice of the torus acting on `G`) into `M` (the weight lattice of the possibly bigger torus acting on the total space of the line bundle).
- `weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}`: A vector containing the weight of the fibre of the line bundle over each vertex of `G`.

# Examples

The trivial line bundle $L\rightarrow \mathbb{P}^2$ with scaling action on each fibre.
Here $T=(\mathbb{C}^\times)^3$ acts on $\mathbb{P}^2$ and $T\times\mathbb{C}^\times$ acts on the total space of $L$, where the 
extra factor $\mathbb{C}^\times$ scales each fiber and preserves the base.
```jldoctest line_bundle
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> GMtoM = ModuleHomomorphism(G.M, M, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]])
GKM vector bundle of rank 1 over GKM graph with 3 nodes and valency 2 with weights:
1: (0, 0, 0, 1)
2: (0, 0, 0, 1)
3: (0, 0, 0, 1)
```
Here is another line bundle on $\mathbb{P}^2$ with a more interesting action than fibrewise scaling:
```jldoctest line_bundle
julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]])
GKM vector bundle of rank 1 over GKM graph with 3 nodes and valency 2 with weights:
1: (1, 0, 0, 0)
2: (0, 1, 0, 0)
3: (0, 0, 1, 0)
```

"""
function line_bundle(
  G::AbstractGKM_graph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Vector{AbstractAlgebra.Generic.FreeModuleElem{R}}
)::GKM_vector_bundle where R<:GKM_weight_type
  return vector_bundle(G, M, GMtoM, reshape(weights, length(weights), 1))
end


@doc raw"""
    vector_bundle(G::AbstractGKM_graph, M::AbstractAlgebra.Generic.FreeModule{R}, GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}, weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}; calculateConnection::Bool=true) -> GKM_vector_bundle

Construct the equivariant vector bundle given by the following datum:#
# Arguments
- `g::G::AbstractGKM_graph`: The GKM graph of the base.
- `M::AbstractAlgebra.Generic.FreeModule{R}`: The weight lattice of the torus acting on the vector bundle.
    This is often bigger than the torus acting on `G`, for example when there is an extra scaling-action on the fibres.
- `GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R}`: The inclusion of `G.M` (the weight lattice of the torus acting on `G`) into `M` (the weight lattice of the possibly bigger torus acting on the total space of the vector bundle).
- `weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}}`: Over each fixed point (i.e., vertex of `G`), the vector bundle splits into a direct sum of $r$ equivariant line bundles, where $r$ is the rank of the vector bundle.
    This argument is a matrix such that `weights[i, j]` is the $j$-th weight at the $i$-th vertex.

# Examples
Let us construct manually (without using `direct_sum()`) the direct sum of the two examples from `line_bundle()`.
```jldoctest
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]]);

julia> V = vector_bundle(G, M, GMtoM, [g[1] g[4]; g[2] g[4]; g[3] g[4]])
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (1, 0, 0, 0), (0, 0, 0, 1)
2: (0, 1, 0, 0), (0, 0, 0, 1)
3: (0, 0, 1, 0), (0, 0, 0, 1)
```
"""
function vector_bundle(
  G::AbstractGKM_graph,
  M::AbstractAlgebra.Generic.FreeModule{R},
  GMtoM::AbstractAlgebra.Generic.ModuleHomomorphism{R},
  weights::Matrix{AbstractAlgebra.Generic.FreeModuleElem{R}};
  calculateConnection::Bool=true
)::GKM_vector_bundle where R<:GKM_weight_type

  @req domain(GMtoM) == G.M "GMtoM must go from gkm.M into M"
  @req rank(kernel(GMtoM)[1]) == 0 "GMtoM must be injective"
  @req codomain(GMtoM) == M "GMtoM must go from gkm.M into M"

  s = size(weights)
  nv = n_vertices(G.g)
  @req s[1] == nv "Weight matrix has wrong dimensions."
  
  for w in weights
    @req parent(w) == M "Weights need to live in M."
  end

  res =  GKM_vector_bundle(G, M, GMtoM, weights, nothing)
  # build connection if it is unique.
  if calculateConnection
    get_connection(res)
  end
  return res
end

@doc raw"""
    rank(V::GKM_vector_bundle) -> Int64

Return the rank of the given GKM vector bundle.

# Example
```jldoctest rank_bdles
julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5));

julia> M = free_module(ZZ, 5);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3], g[4]]);

julia> L = line_bundle(G, M, GMtoM, [g[5], g[5], g[5], g[5]]);

julia> rank(direct_sum(L, L, L))
3
```

"""
function Oscar.rank(V::GKM_vector_bundle)::Int64
  return size(V.w)[2]
end

@doc raw"""
    tangent_bd(G::AbstractGKM_graph; scaling_weight::Int64 = 1) -> GKM_vector_bundle

Return the tangent bundle of `G`. The torus is enlarged by one dimension where the extra factor scales the fibers of the tangent bundle.
The default weight is 1, but can be changed using the optional argument `scaling_weight` if desired.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> T = tangent_bd(G)
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (1, -1, 0, 1), (1, 0, -1, 1)
2: (-1, 1, 0, 1), (0, 1, -1, 1)
3: (-1, 0, 1, 1), (0, -1, 1, 1)

julia> P = projectivization(T)
GKM graph is valid but not 3-independent, so connections may not be unique.
GKM graph with 6 nodes, valency 3 and axial function:
[1]_2 -> [1]_1 => (0, -1, 1, 0)
[2]_1 -> [1]_1 => (-1, 1, 0, 0)
[2]_2 -> [1]_2 => (-1, 1, 0, 0)
[2]_2 -> [2]_1 => (-1, 0, 1, 0)
[3]_1 -> [1]_2 => (-1, 0, 1, 0)
[3]_1 -> [2]_1 => (0, -1, 1, 0)
[3]_2 -> [1]_1 => (-1, 0, 1, 0)
[3]_2 -> [2]_2 => (0, -1, 1, 0)
[3]_2 -> [3]_1 => (-1, 1, 0, 0)

julia> betti_numbers(P)
4-element Vector{Int64}:
 1
 2
 2
 1
```
"""
function tangent_bd(G::AbstractGKM_graph; scaling_weight::Int64 = 1)::GKM_vector_bundle

  return _co_tangent_bundle(G, scaling_weight, 1)
end

@doc raw"""
    cotangent_bd(G::AbstractGKM_graph; scaling_weight::Int64 = 1) -> GKM_vector_bundle

Return the cotangent bundle of `G`. The torus is enlarged by one dimension where the extra factor scales the fibers of the tangent bundle.
The default weight is 1, but can be changed using the optional argument `scaling_weight` if desired.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 3)
GKM graph with 4 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
4 -> 1 => (-1, 0, 0, 1)
4 -> 2 => (0, -1, 0, 1)
4 -> 3 => (0, 0, -1, 1)

julia> T = cotangent_bd(G)
GKM vector bundle of rank 3 over GKM graph with 4 nodes and valency 3 with weights:
1: (-1, 1, 0, 0, 1), (-1, 0, 1, 0, 1), (-1, 0, 0, 1, 1)
2: (1, -1, 0, 0, 1), (0, -1, 1, 0, 1), (0, -1, 0, 1, 1)
3: (1, 0, -1, 0, 1), (0, 1, -1, 0, 1), (0, 0, -1, 1, 1)
4: (1, 0, 0, -1, 1), (0, 1, 0, -1, 1), (0, 0, 1, -1, 1)
```
"""
function cotangent_bd(G::AbstractGKM_graph; scaling_weight::Int64 = 1)::GKM_vector_bundle

  return _co_tangent_bundle(G, scaling_weight, -1)
end

function _co_tangent_bundle(G::AbstractGKM_graph, scaling_weight::Int64, duality::Int64)::GKM_vector_bundle
  R = base_ring(G.M)
  nv = n_vertices(G.g)
  r = rank_torus(G)
  M = free_module(R, r+1)
  g = gens(M)
  GMtoM = ModuleHomomorphism(G.M, M, [g[i] for i in 1:r])
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(R))}}(undef, nv, valency(G))
  for v in 1:nv
    ct = 0
    for w in all_neighbors(G.g, v)
      ct += 1
      weightMatrix[v, ct] = duality * GMtoM(G.w[Edge(v, w)]) + scaling_weight * g[r+1]
    end
  end
  return vector_bundle(G, M, GMtoM, weightMatrix)
end

@doc raw"""
    get_connection(V::GKM_vector_bundle)

Return the connection of the given vector bundle, if it is unique or has been set manually.
If the vector bundle does not admit a unique connection and it has not bene set manually, return `nothing`.

# Mathematical description:
This is the same concept as a [Connection](Connections.md) on a GKM graph.
Let `G` be the GKM graph that is the basis of the vector bundle `V`.
Assume that `G` comes from a GKM variety $X$.
Then each edge of `e` corresponds to an invariant rational curve $C_e$ in $X$.
If $V$ is an equivariant line bundle overe $X$, then its restriction to $C_e\cong\mathbb{P}^1$ splits into a direct sum of equivariant line bundles.
This defines a bijection between the direct summands of $V$ at $\text{src}(e)$ and the direct summands of $V$ at $\text{dst}(e)$.

This bijection is recorded in the returned object.

# Example
```julia-repl
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]]);

julia> V = vector_bundle(G, M, GMtoM, [g[1] g[4]; g[2] g[4]; g[3] g[4]])
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (1, 0, 0, 0), (0, 0, 0, 1)
2: (0, 1, 0, 0), (0, 0, 0, 1)
3: (0, 0, 1, 0), (0, 0, 0, 1)

julia> get_connection(V)
Dict{Tuple{Edge, Int64}, Int64} with 12 entries:
  (Edge(3, 1), 2) => 2
  (Edge(1, 2), 1) => 1
  (Edge(3, 1), 1) => 1
  (Edge(1, 2), 2) => 2
  (Edge(3, 2), 1) => 1
  (Edge(2, 3), 1) => 1
  (Edge(3, 2), 2) => 2
  (Edge(2, 3), 2) => 2
  (Edge(1, 3), 1) => 1
  (Edge(2, 1), 1) => 1
  (Edge(1, 3), 2) => 2
  (Edge(2, 1), 2) => 2
```
It is visible here that the vector bundle is a direct sum of two line bundles, since we have `(e, i) => i` for each edge `e` and index `i`.
The output will be more complicated when the vector bundle does not split into line bundles.
"""
function get_connection(V::GKM_vector_bundle)
  
  if isnothing(V.con)
    V.con = _build_vector_bundle_connection(V)
  end
  return V.con
end

# Return the unique GKM conncetion of the vector bundle or nothing if it is not uniquely determined.
function _build_vector_bundle_connection(V::GKM_vector_bundle)

  con = Dict{Tuple{Edge, Int64}, Int64}()
  weights = V.w

  G = V.gkm
  rk = rank(V)

  for e in edges(G.g)
    @req !is_zero(G.w[e]) "Weight zero edge found."
    v = src(e)
    w = dst(e)
    we = V.GMtoM(G.w[e])
    for i in 1:rk
      wi = weights[v, i]
      haveFoundJ = false
      for j in 1:rk
        wj = weights[w, j]
        wdif = wi - wj
        if rank(matrix([ wdif; we ])) == 1 # if true, (v,i) belongs to (w,j)
          if haveFoundJ
            # connection is not unique, so return nothing.
            return nothing
          else
            # have found a unique (so far) candidate for j.
            con[(e, i)] = j
            con[(reverse(e), j)] = i
            haveFoundJ = true
          end
        end
      end
      if !haveFoundJ
        return nothing
      end
    end
  end
  return con
end

@doc raw"""
    direct_sum(V::GKM_vector_bundle{R}...) -> GKM_vector_bundle

Return the direct sum of the given vector bundles.
This requires all bundles to have the same base GKM graph and the same character lattice.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]]);

julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = direct_sum(V1, V2)
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (0, 0, 0, 1), (1, 0, 0, 0)
2: (0, 0, 0, 1), (0, 1, 0, 0)
3: (0, 0, 0, 1), (0, 0, 1, 0)
```
"""
function direct_sum(V::GKM_vector_bundle{R}...)::GKM_vector_bundle where R<:GKM_weight_type
  n = length(V)
  @req n >= 1 "Need at least one direct summand."
  for i in 1:n, j in 1:n
    @req V[i].gkm == V[j].gkm "Vector bundles need to have the same GKM base."
    @req V[i].M == V[j].M "Vector bundles need to have the same character lattice."
    @req V[i].GMtoM == V[j].GMtoM "V.GMtoM needs to be constant among direct summands."
  end
  G = V[1].gkm
  M = V[1].M
  GMtoM = V[1].GMtoM
  weights = hcat((V[i].w for i in 1:n)...)
  res = vector_bundle(G, M, GMtoM, weights)

  # infer connection from direct summands
  if isnothing(get_connection(res)) && !any([isnothing(get_connection(V[i])) for i in 1:n])
    con = Dict{Tuple{Edge, Int64}, Int64}()
    offset = 0
    for i in 1:n
      conI = get_connection(V[i])
      for k in keys(conI)
        e = k[1]
        a = k[2]
        b = conI[k]
        con[(e, a+offset)] = b+offset
        con[(reverse(e), b+offset)] = a+offset
      end
      offset += rank(V[i])
    end
    res.con = con
  end
  return res
end


function Base.show(io::IO, V::GKM_vector_bundle)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "GKM vector bundle")
  else
    # nested printing allowed, preferably terse
    print(io, "GKM vector bundle of rank $(rank(V)) over GKM graph with $(n_vertices(V.gkm.g)) vertices")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", V::GKM_vector_bundle)

  print(io, "GKM vector bundle of rank $(rank(V)) over $(V.gkm) with weights:")
  rk = rank(V)
  for v in 1:n_vertices(V.gkm.g)
    print(io, "\n$(V.gkm.labels[v]): ")
    for i in 1:rk
      print(io, V.w[v,i])
      if i<rk
        print(io, ", ")
      end
    end
  end
end

@doc raw"""
    dual(V::GKM_vector_bundle) -> GKM_vector_bundle

Return the dual equivariant vector bundle.

# Example
```jldoctest dual_vector_bundles
julia> G = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, 5));

julia> M = free_module(ZZ, 5);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3], g[4]]);

julia> L = line_bundle(G, M, GMtoM, [g[5], g[5], g[5], g[5]])
GKM vector bundle of rank 1 over GKM graph with 4 nodes and valency 2 with weights:
1: (0, 0, 0, 0, 1)
2: (0, 0, 0, 0, 1)
3: (0, 0, 0, 0, 1)
4: (0, 0, 0, 0, 1)

julia> dual(L)
GKM vector bundle of rank 1 over GKM graph with 4 nodes and valency 2 with weights:
1: (0, 0, 0, 0, -1)
2: (0, 0, 0, 0, -1)
3: (0, 0, 0, 0, -1)
4: (0, 0, 0, 0, -1)

```
"""
function Oscar.dual(V::GKM_vector_bundle)::GKM_vector_bundle
  res = vector_bundle(V.gkm, V.M, V.GMtoM, -V.w; calculateConnection=false)
  res.con = V.con
  return res
end

@doc raw"""
    projectivization(V::GKM_vector_bundle) -> AbstractGKM_graph

Return the projectivisation of the given equivariant vector bundle.
!!! note
    If the given bundle does not admit a unique connection, it must be specified manually by setting the field `V.con`.

# Example
```jldoctest projectivization
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]]);

julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = direct_sum(V1, V2)
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (0, 0, 0, 1), (1, 0, 0, 0)
2: (0, 0, 0, 1), (0, 1, 0, 0)
3: (0, 0, 0, 1), (0, 0, 1, 0)

julia> P = projectivization(V)
GKM graph with 6 nodes, valency 3 and axial function:
[1]_2 -> [1]_1 => (-1, 0, 0, 1)
[2]_1 -> [1]_1 => (-1, 1, 0, 0)
[2]_2 -> [1]_2 => (-1, 1, 0, 0)
[2]_2 -> [2]_1 => (0, -1, 0, 1)
[3]_1 -> [1]_1 => (-1, 0, 1, 0)
[3]_1 -> [2]_1 => (0, -1, 1, 0)
[3]_2 -> [1]_2 => (-1, 0, 1, 0)
[3]_2 -> [2]_2 => (0, -1, 1, 0)
[3]_2 -> [3]_1 => (0, 0, -1, 1)
```
The naming convention for the vertices of the projectivization's GKM graph is `[v]_i` where `v` is a vertex of the original GKM graph and `i` is the index
of the line bundle direct summand to which this vertex of the projectivization corresponds.
"""
function Oscar.projectivization(V::GKM_vector_bundle)::AbstractGKM_graph
  con = get_connection(V)
  @req !isnothing(con) "GKM vector bundle needs connection for projectivization."
  G = V.gkm
  nv = n_vertices(G.g)
  rk = rank(V)

  Gres = Graph{Undirected}(nv * rk)
  # labels = String[]
  # sizehint!(labels, nv * rk)
  # #build labels
  # for v in 1:nv
  #   for i in 1:rk
  #     push!(labels, "[" * G.labels[v] * "]_$i")
  #   end
  # end
  labels = Vector{String}(undef, nv * rk)
  for v in 1:nv
    for i in 1:rk
      labels[(v-1)*rk + i] = "[$(G.labels[v])]_$i"
    end
  end


  weightType = typeof(_get_weight_type(G))
  res = gkm_graph(Gres, labels, V.M, Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{weightType}}(); checkLabels=false)

  Gcon = get_connection(G)
  resCon = Dict{Tuple{Edge, Edge}, Edge}()

  # add edges corresponding to original edges
  for e in edges(G.g)
    v = src(e)
    w = dst(e)
    for i in 1:rk
      vInd = (v-1)*rk + i
      j = con[(e, i)]
      wInd = (w-1)*rk + j
      add_edge!(res, vInd, wInd, V.GMtoM(G.w[e]))

      if !isnothing(Gcon)
        eNew = Edge(vInd, wInd)
        for u in all_neighbors(G.g, v)
          uInd = (u-1)*rk + con[(Edge(v, u), i)]
          ei = Edge(vInd, uInd)
          up = dst(Gcon.con[(e, Edge(v, u))])
          upInd = (up - 1)*rk + con[(Edge(w, up), j)]
          epi = Edge(wInd, upInd)
          resCon[(eNew, ei)] = epi
          resCon[(reverse(eNew), epi)] = ei
        end
        for k in 1:rk
          i == k && continue
          kInd = (v-1)*rk + k
          l = con[(e, k)]
          lInd = (w-1)*rk + l
          ei = Edge(vInd, kInd)
          epi = Edge(wInd, lInd)
          resCon[(eNew, ei)] = epi
          resCon[(reverse(eNew), epi)] = ei
        end
      end
    end
  end
  # add edges over each vertex
  for v in 1:nv
    for i in 1:rk
      for j in (i+1):rk
        viInd = (v-1)*rk + i
        vjInd = (v-1)*rk + j
        wNew = V.w[v, j] - V.w[v, i]
        @req !iszero(wNew) "Vector bundle has two identical weights over vertex $v (indices $i, $j)"
        add_edge!(res, viInd, vjInd, wNew)

        if !isnothing(Gcon)
          e = Edge(viInd, vjInd)
          resCon[(e, e)] = reverse(e)
          resCon[(reverse(e), reverse(e))] = e
          for k in 1:rk
            (i == k || j == k) && continue
            vkInd = (v-1)*rk + k
            ei = Edge(viInd, vkInd)
            epi = Edge(vjInd, vkInd)
            resCon[(e, ei)] = epi
            resCon[(reverse(e), epi)] = ei
          end
          for w in all_neighbors(G.g, v)
            eDown = Edge(v, w)
            k = con[(eDown, i)]
            l = con[(eDown, j)]
            wkInd = (w-1)*rk + k
            wlInd = (w-1)*rk + l
            ei = Edge(viInd, wkInd)
            epi = Edge(vjInd, wlInd)
            resCon[(e, ei)] = epi
            resCon[(reverse(e), epi)] = ei
          end
        end
      end
    end
  end

  if !isnothing(Gcon)
    resConObj = build_GKM_connection(res, resCon)
    set_connection!(res, resConObj)
  end

  if !isvalid(res)
    println("Warning: resulting projective bundle is not a valid GKM graph (see reason above).")
  end

  return res
end

function _calculate_connection_a(V::GKM_vector_bundle; check::Bool=true)
  has_attribute(V, :connectionA) && return

  rV = rank(V)
  connectionA = Dict{Tuple{Edge, Int64}, ZZRingElem}()
  con = get_connection(V)
  @req !isnothing(con) "V needs a connection to calculate connection a's!"

  for e in edges(V.gkm.g)
    eW = V.gkm.w[e]
    for i in 1:rV
      wei = V.w[src(e), i]
      k = con[(e, i)]
      wepi = V.w[dst(e), k]
      wdif = wei - wepi

      if check
        @req rank(matrix([ wdif; eW ])) == 1 "connection of vector bundle is incompatible with GKM graph"
      end

      ai::ZZRingElem = ZZ(0)

      for j in 1:rank(V.gkm.M)
        if eW[j] != 0
          tmp = wdif[j] // eW[j]
          @req denominator(tmp) == 1 "GKM connection's a_i's must be integers!" # Assumption: x//y is integer if and only if denominator(x//y) == 1 in Oscar.
          ai = ZZ(tmp)
          break
        end
      end

      connectionA[(e, i)] = ai
      connectionA[(reverse(e), k)] = ai
    end
  end
  set_attribute!(V, :connectionA, connectionA)
end

function _calculate_weight_classes(V::GKM_vector_bundle)
  G = V.gkm
  has_attribute(V, :normalClasses) && has_attribute(V, :weightClasses) && return

  # This here needs revision later.
  @req G.M == V.M "Weight classes are currently only supported for G.M == V.M and GMtoM = identity."

  nv = n_vertices(G.g)
  rV = rank(V)
  rT = rank_torus(G)
  R = G.equivariantCohomology.coeffRing
  t = gens(R)
  weightClasses = Matrix{QQMPolyRingElem}(undef, nv, rV)
  normalClasses = Vector{QQMPolyRingElem}(undef, nv)
  
  # normal and weight classes:
  for v in 1:nv
    nc = one(R)
    for i in 1:rV
      w = zero(R)
      for j in 1:rT
        w += V.w[v, i][j] * t[j]
      end
      weightClasses[v, i] = w
      nc = nc * w
    end
    normalClasses[v] = nc
  end

  set_attribute!(V, :normalClasses, normalClasses)
  set_attribute!(V, :weightClasses, weightClasses)
end

# This only works if _calculate_weight_classes(V) was called before!
function _fiber_normal_weight(v::Int64, V::GKM_vector_bundle)
  return get_attribute(V, :normalClasses)[v]
end

# This only works if _calculate_weight_classes(V) was called before!
function _fiber_summand_weight(v::Int64, i::Int64, V::GKM_vector_bundle)
  return get_attribute(V, :weightClasses)[v,i]
end

# This only works if _calculate_connection_a(V) was called before!
function _fiber_connection_a(e::Edge, i::Int64, V::GKM_vector_bundle)
  return get_attribute(V, :connectionA)[(e, i)]
end