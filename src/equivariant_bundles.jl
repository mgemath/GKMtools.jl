export wedge_product, sym_product, baseof, gkm_line_bundle_of_toric, gkm_vector_bundle_of_toric

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

Construct the equivariant vector bundle given by the following datum:
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
    baseof(V::GKM_vector_bundle) -> AbstractGKM_graph

Return the base of the given GKM vector bundle.
"""
function baseof(V::GKM_vector_bundle)::AbstractGKM_graph
  return V.gkm
end

@doc raw"""
    total_space(V::GKM_vector_bundle) -> AbstractGKM_graph

Return the total space of the given GKM vector bundle as a GKM graph.

The total space is constructed by adding `rank(V)` standalone flags at each vertex of the base,
with weights given by the fibre weights transformed via `GMtoM`.

If the base has a connection set, and the vector bundle has a connection, the total space
connection is also constructed. Similarly, if `H2` is computed for the base (and the optional argument
`copy_curve_classes` was not manually set to `false`), it is copied
to the total space.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 2)
GKM graph with 3 nodes, valency 2 and axial function:
2 -> 1 => (-1, 1, 0)
3 -> 1 => (-1, 0, 1)
3 -> 2 => (0, -1, 1)

julia> M = free_module(ZZ, 4);

julia> GMtoM = ModuleHomomorphism(G.M, M, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]])
GKM vector bundle of rank 1 over GKM graph with 3 nodes and valency 2 with weights:
1: (0, 0, 0, 1)
2: (0, 0, 0, 1)
3: (0, 0, 0, 1)

julia> T = total_space(V)
GKM graph with 3 nodes, valency 3 and axial function:
2 -> 1 => (-1, 1, 0, 0)
3 -> 1 => (-1, 0, 1, 0)
3 -> 2 => (0, -1, 1, 0)
Standalone flags:
1.3 => (0, 0, 0, 1)
2.3 => (0, 0, 0, 1)
3.3 => (0, 0, 0, 1)

julia> valency(T)
3

julia> is_compact(T)
false
```

Another well-known example is the total space of the cotangent bundle, which is Calabi--Yau.
```jldoctest
julia> F = gkm_3d_twisted_flag()
GKM graph with 6 nodes, valency 3 and axial function:
2 -> 1 => (0, -1)
3 -> 2 => (1, 0)
4 -> 1 => (1, -2)
4 -> 3 => (-1, 1)
5 -> 2 => (1, -1)
5 -> 4 => (0, -1)
6 -> 1 => (1, -1)
6 -> 3 => (2, -1)
6 -> 5 => (1, 0)

julia> X = total_space(cotangent_bd(F));

julia> print_curve_classes(X)
2 -> 1: (0, 1), Chern number: 0
3 -> 2: (-1, 1), Chern number: 0
4 -> 1: (1, 0), Chern number: 0
4 -> 3: (-2, 1), Chern number: 0
5 -> 2: (1, 0), Chern number: 0
5 -> 4: (-1, 1), Chern number: 0
6 -> 1: (1, 1), Chern number: 0
6 -> 3: (1, 0), Chern number: 0
6 -> 5: (0, 1), Chern number: 0

```
"""
function Oscar.total_space(V::GKM_vector_bundle{R}; copy_curve_classes::Bool=true)::AbstractGKM_graph{R} where R <: GKM_weight_type
  base = V.gkm
  r = rank(V)
  nv = n_vertices(base.g)
  base_val = valency(base)

  # Build connection if both base and bundle have connections
  base_con = get_connection(base)
  bundle_con = get_connection(V)  # This is V.con: Dict{Tuple{Edge, Int64}, Int64}

  # start with a deep copy of the base, with substitited weights according to V.GMtoM
  # This also copies the connection.
  total = substitute_torus(base, V.GMtoM)

  # Add standalone flags for fibres at each vertex
  for v in 1:nv
    for i in 1:r
      add_standalone_flag!(total, v, V.w[v, i])
    end
  end

  # Copy H2 if it exists for the base
  if copy_curve_classes && !isnothing(base.curveClasses)
    # create copy of H2 object and update its reference to parent GKM graph and Chern numbers.
    H2_copy = deepcopy(base.curveClasses)
    H2 = H2_copy.H2
    edgeLattice = H2_copy.edgeLattice
    quotientMap = H2_copy.quotientMap
    edgeToGenIndex = H2_copy.edgeToGenIndex
    
    # chern numbers may have changed
    dualConeRaySum, C, H2ToCN = _finish_GKM_H2(edgeLattice, H2, quotientMap, total, edgeToGenIndex)
    newH2 = GKM_H2(total, edgeLattice, H2, edgeToGenIndex, quotientMap, dualConeRaySum, C, H2ToCN, nothing, nothing) # the last two nothings forget Seidel space data.
    total.curveClasses = newH2
  end

  if !isnothing(base_con) && !isnothing(bundle_con)
    # This total connection is already correct on the base.
    @req !isnothing(total.connection) "substitute_torus failed to copy connection from base space."
    total_con = total.connection.con

    for e in keys(total_con)
      append!(total_con[e], zeros(Int64, r))
    end

    # Add fiber data to connection.
    for (e, i) in keys(bundle_con)
      total_con[e][base_val + i] = base_val + bundle_con[(e, i)]
    end

    # Use build_GKM_connection to create the full connection (this computes a-values)
    total_con_obj = build_GKM_connection(total, total_con)
    set_connection!(total, total_con_obj)
  else
    # Make sure we don't keep a half-defined connection from the base installed.
    # Whithout this, we could get the connection part from the base but nothing on the fibers of V.
    total.connection = nothing
  end

  return total
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
  val = valency(G)
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(zero(R))}}(undef, nv, val)

  # For each vertex, iterate over all flags (not just edge neighbors)
  for v in 1:nv
    for i in 1:val
      # Use the flag-based weight structure
      weightMatrix[v, i] = duality * GMtoM(G.weights_at_vertex[v][i]) + scaling_weight * g[r+1]
    end
  end

  # Create the vector bundle
  V = vector_bundle(G, M, GMtoM, weightMatrix; calculateConnection = false)

  # If the base graph has a connection, compute the induced connection on the tangent/cotangent bundle
  base_con = get_connection(G)
  if !isnothing(base_con)
    # The tangent bundle's fiber flags correspond directly to the base graph's flags
    # So the connection on the tangent bundle is induced from the base connection
    bundle_con = Dict{Tuple{Edge, Int64}, Int64}()

    for e in edges(G.g)
      for i in 1:val
        # The i-th fiber flag at src(e) connects to con[e][i]-th fiber flag at dst(e)
        j = base_con.con[e][i]
        bundle_con[(e, i)] = j
        bundle_con[(reverse(e), j)] = i
      end
    end

    V.con = bundle_con
    @req isvalid(V.con, V) "Co/Tangent bundle induces invalid connection on V."
  end

  return V
end

@doc raw"""
    get_connection(V::GKM_vector_bundle)

Return the connection of the given vector bundle, if it is unique or has been set manually.
If the vector bundle does not admit a unique connection and it has not bene set manually, return `nothing`.

# Mathematical description:
This is the same concept as a [Connection](Connections.md) on a GKM graph.
Let `G` be the GKM graph that is the basis of the vector bundle `V`.
Assume that `G` comes from a GKM space $X$.
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

@doc raw"""
    get_any_connection(V::GKM_vector_bundle)

Return any compatible connection for the given vector bundle, if one exists.
The result does not necessarily agree with the splittings of $V$ over the $T$-stable $\mathbb{P}^1$s.
However, for many applications in Gromov--Witten theory, any compatible connection is enough.
"""
function get_any_connection(V::GKM_vector_bundle)
  con = get_connection(V)
  if !isnothing(con)
    return con
  elseif isnothing(V.anyConnection)
    V.anyConnection = _build_any_vector_bundle_connection(V)
  end
  return V.anyConnection
end

@doc raw"""
    isvalid(con::Dict{Tuple{Edge, Int64}, Int64}, V::GKM_vector_bundle; printDiagnostics::Bool=true) -> Bool

Check if a vector bundle connection is valid.

A connection is valid if:
1. It has entries for all edges (in both directions)
2. The connection respects the involution: if `(e,i)` maps to `j`, then `(reverse(e),j)` maps to `i`
3. For each edge `e` and fiber index `i`, there exists an integer `a_i` such that:
   `w[dst(e), con[(e,i)]] = w[src(e), i] - a_i * edge_weight(e)`

# Arguments
- `con`: The connection dictionary mapping (edge, fiber_index) to connected fiber index
- `V`: The vector bundle
- `printDiagnostics`: If true, print diagnostic messages when validation fails

# Example
```julia-repl
julia> G = projective_space(GKM_graph, 2);

julia> TG = tangent_bd(G);

julia> con = get_any_connection(TG);

julia> isvalid(con, TG)
true
```
"""
function isvalid(con::Dict{Tuple{Edge, Int64}, Int64}, V::GKM_vector_bundle; printDiagnostics::Bool=true)::Bool

  G = V.gkm
  rk = rank(V)

  # Check that all edges have connection entries
  for e in edges(G.g)
    # Check forward direction
    for i in 1:rk
      if !haskey(con, (e, i))
        printDiagnostics && println("Connection missing entry for edge $e, fiber $i")
        return false
      end

      j = con[(e, i)]

      # Check that j is in valid range
      if j < 1 || j > rk
        printDiagnostics && println("Connection maps edge $e, fiber $i to invalid fiber index $j (rank is $rk)")
        return false
      end
    end

    # Check reverse direction
    for j in 1:rk
      if !haskey(con, (reverse(e), j))
        printDiagnostics && println("Connection missing entry for edge $(reverse(e)), fiber $j")
        return false
      end
    end
  end

  # Check that connection respects involution
  for e in edges(G.g)
    for i in 1:rk
      j = con[(e, i)]
      i_back = con[(reverse(e), j)]
      if i_back != i
        printDiagnostics && println("Connection involution violated: con[$e][$i] = $j, but con[$(reverse(e))][$j] = $i_back ≠ $i")
        return false
      end
    end
  end

  # Check that a-values exist and are integers
  # We need to check both e and reverse(e) to ensure consistency
  for e_base in edges(G.g)
    for edge in [e_base, reverse(e_base)]
      eW = V.GMtoM(_w(G, edge))

      for i in 1:rk
        j = con[(edge, i)]
        wi = V.w[src(edge), i]
        wj = V.w[dst(edge), j]
        wdif = wi - wj

        # Check that wdif and eW are linearly dependent (rank = 1 or both zero)
        mat_rank = rank(matrix([wdif; eW]))
        if mat_rank > 1
          printDiagnostics && println("Connection incompatible with weights at edge $edge, fiber $i: w_diff and edge_weight are not linearly dependent")
          return false
        end

        # Find a_i such that wdif = a_i * eW
        ai::Union{Nothing, ZZRingElem} = nothing

        for k in 1:rank(V.M)
          if eW[k] != 0
            tmp = wdif[k] // eW[k]
            if denominator(tmp) != 1
              printDiagnostics && println("Connection a-value is not an integer at edge $edge, fiber $i: wdif[$k]/eW[$k] = $(wdif[k])/$(eW[k]) = $tmp")
              return false
            end
            ai = ZZ(tmp)
            break
          end
        end

        # If eW is zero, wdif must also be zero
        if isnothing(ai)
          if !iszero(wdif)
            printDiagnostics && println("Edge weight is zero at edge $edge, but weight difference is non-zero for fiber $i")
            return false
          end
          # If both are zero, any a_i would work, so we consider this valid
        else
          # Verify that a_i * eW = wdif for ALL components
          if G.weightType(ai) * eW != wdif
            printDiagnostics && println("Connection a-value inconsistent at edge $edge, fiber $i: $ai * eW ≠ wdif")
            return false
          end
        end
      end
    end
  end

  return true
end

# Return the unique GKM conncetion of the vector bundle or nothing if it is not uniquely determined.
function _build_vector_bundle_connection(V::GKM_vector_bundle)

  con = Dict{Tuple{Edge, Int64}, Int64}()
  weights = V.w

  G = V.gkm
  rk = rank(V)

  for e in edges(G.g)
    @req !is_zero(_w(G, e)) "Weight zero edge found."
    v = src(e)
    w = dst(e)
    we = V.GMtoM(_w(G, e))
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

function _build_any_vector_bundle_connection(V::GKM_vector_bundle)

  con = Dict{Tuple{Edge, Int64}, Int64}()
  weights = V.w

  G = V.gkm
  rk = rank(V)

  for e in edges(G.g)
    #println("e=$e")
    @req !is_zero(_w(G, e)) "Weight zero edge found."
    v = src(e)
    w = dst(e)
    we = V.GMtoM(_w(G, e))
    #println("we = $we")
    # make sure not to allocate some epi to more than one ei.
    allocatedJs = Vector{Int64}()
    for i in 1:rk
      #println("  i=$i")
      wi = weights[v, i]
      for j in 1:rk
        #println("    j=$j")
        wj = weights[w, j]
        wdif = wi - wj
        #println("    wdif=$wdif")
        if rank(matrix([ wdif; we ])) == 1 && !(j in allocatedJs)
          
          # j is only a candidate for i if the resulting ai is an integer.
          aiIntegral::Bool = false
          for k in 1:rank(V.M)
            if we[k] != 0
              tmp = wdif[k] // we[k]
              aiIntegral = denominator(tmp) == 1
              #println("    aiIntegral = $aiIntegral")
              break
            end
          end
          !aiIntegral && continue

          # have found a match for i.
          #println("  set j = $j")
          con[(e, i)] = j
          con[(reverse(e), j)] = i
          push!(allocatedJs, j)
          break
        end
      end
      if !haskey(con, (e, i))
        println("No connection image found for ($e, $i)! The GKM vector bundle does not admit a connection.")
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
  for i in 2:n
    if !(V[i].gkm == V[1].gkm)
      @warn "Vector bundles could be defined on different GKM bases."
    end
  # @req V[i].gkm == V[j].gkm "Vector bundles need to have the same GKM base."
  @req V[i].M == V[1].M "Vector bundles need to have the same character lattice."
  @req V[i].GMtoM == V[1].GMtoM "V.GMtoM needs to be constant among direct summands."
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

@doc raw"""
    +(V::GKM_vector_bundle, W::GKM_vector_bundle) -> GKM_vector_bundle

Return the direct sum of two vector bundles.
This requires all bundles to have the same base GKM graph and the same character lattice.

# Example
```jldoctest
julia> G = projective_space(GKM_graph, 2);

julia> M = free_module(ZZ, 4);

julia> g = gens(M);

julia> GMtoM = ModuleHomomorphism(G.M, M, [g[1], g[2], g[3]]);

julia> V1 = line_bundle(G, M, GMtoM, [gens(M)[4], gens(M)[4], gens(M)[4]]);

julia> V2 = line_bundle(G, M, GMtoM, [gens(M)[1], gens(M)[2], gens(M)[3]]);

julia> V = V1 + V2
GKM vector bundle of rank 2 over GKM graph with 3 nodes and valency 2 with weights:
1: (0, 0, 0, 1), (1, 0, 0, 0)
2: (0, 0, 0, 1), (0, 1, 0, 0)
3: (0, 0, 0, 1), (0, 0, 1, 0)
```
"""
function +(V::GKM_vector_bundle, W::GKM_vector_bundle)::GKM_vector_bundle
  return direct_sum(V, W)
end

@doc raw"""
    det(V::GKM_vector_bundle) -> GKM_vector_bundle

Return the determinant of the vector bundle `V`, that is $\wedge^{\mathrm{rank}(V)} V$.
"""
function det(V::GKM_vector_bundle)::GKM_vector_bundle
  return wedge_product(V, rank(V))
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

# Examples
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

!!! warning
    This function may create
    GKM graphs that fail `isvalid` because they are not 2-independent.
    This happens in the following example of $\mathbb{P}(T_X\oplus T^*_X)$ with $X=\mathbb{P}^1$,
    because the differences of fiber weights in the tangent and cotangent bundle
    coincide up to sign.

```jldoctest
julia> P1 = projective_space(GKM_graph, 1)
GKM graph with 2 nodes, valency 1 and axial function:
2 -> 1 => (-1, 1)

julia> T = tangent_bd(P1) + cotangent_bd(P1)
GKM vector bundle of rank 2 over GKM graph with 2 nodes and valency 1 with weights:
1: (1, -1, 1), (-1, 1, 1)
2: (-1, 1, 1), (1, -1, 1)

julia> P = projectivization(T)
GKM graph is not 2-independent.
┌ Warning: Creating GKM subgraph of invalid gkm graph. This may result in undefined behavior.
└ @ GKMtools GKMsubgraphs.jl:243
GKM graph with 4 nodes, valency 2 and axial function:
[1]_2 -> [1]_1 => (2, -2, 0)
[2]_1 -> [1]_1 => (-1, 1, 0)
[2]_2 -> [1]_2 => (-1, 1, 0)
[2]_2 -> [2]_1 => (-2, 2, 0)
```
"""
function Oscar.projectivization(V::GKM_vector_bundle)::AbstractGKM_graph

  # Nota bene: the correctness of this function is depends crucially on the implementation
  # details of total_space and blow_up. Specifically, it assumes that:
  #   1. total_space puts all base flags first and then all fiber flags at each vertex.
  #   2. blowup creates first all exceptional vertices for subgraph vertex 1, then all
  #       for subgraph vertex 2, and so on.

  @req !isnothing(get_connection(V)) "GKM vector bundle needs connection for projectivization."

  base = V.gkm
  nv = n_vertices(base.g)
  rk = rank(V)
  val_base = valency(base)

  # Step 1: Construct the total space
  total = GKMtools.total_space(V)

  # Step 2: Create subgraph of all vertices with all non-fiber flags
  # The base flags are indices 1..val_base (by implementation of total_space)
  # The fiber flags are indices (val_base+1)..(val_base+rk)
  # We want to blow up at the subgraph that includes all vertices but only base flags

  base_flags_at_vertices = Vector{Vector{Int64}}(undef, nv)
  for v in 1:nv
    base_flags_at_vertices[v] = collect(1:val_base)
  end

  # Create subgraph from all vertices with only base flags (not fiber flags)
  base_subgraph = gkm_subgraph_from_flags(total, collect(1:nv), base_flags_at_vertices)

  # Step 3: Blow up at this subgraph
  blowup_result = blow_up(base_subgraph)

  # Step 4: Change blowup labels to projectivization labels
  labels = Vector{String}(undef, nv * rk)
  ctr = 0
  for i in 1:nv
    for j in 1:rk
      ctr += 1
      labels[ctr] = "[$i]_$j"
    end
  end
  blowup_result.self.labels = labels

  # Step 5: Return the exceptional locus (the subgraph part of the blowup)
  return blowup_result.self
end

function _calculate_connection_a(V::GKM_vector_bundle; check::Bool=true)
  has_attribute(V, :connectionA) && return

  rV = rank(V)
  connectionA = Dict{Tuple{Edge, Int64}, ZZRingElem}()
  con = get_any_connection(V)
  @req !isnothing(con) "V needs a connection to calculate connection a's!"

  for e in edges(V.gkm.g)
    eW = V.GMtoM(_w(V.gkm, e))
    for i in 1:rV
      wei = V.w[src(e), i]
      k = con[(e, i)]
      wepi = V.w[dst(e), k]
      wdif = wei - wepi

      if check
        @req rank(matrix([ wdif; eW ])) == 1 "connection of vector bundle is incompatible with GKM graph"
      end

      # Find ai such that wdif = ai * eW
      # We need to find a non-zero component of eW
      ai::Union{Nothing, ZZRingElem} = nothing

      for j in 1:rank(V.gkm.M)
        if eW[j] != 0
          tmp = wdif[j] // eW[j]
          @req denominator(tmp) == 1 "GKM connection's a_i's must be integers!" # Assumption: x//y is integer if and only if denominator(x//y) == 1 in Oscar.
          ai = ZZ(tmp)
          break
        end
      end

      if isnothing(ai)
        # eW is zero, so wdif must also be zero (due to rank check)
        # In this case, ai is not well-defined - the connection a-value is arbitrary
        # This should not happen in a well-formed GKM graph
        error("Edge weight is zero for edge $e - cannot compute connection a-value")
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
  @req G.M === V.M "Weight classes are currently only supported for G.M === V.M and GMtoM = identity."

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
# If the connection is not unique, any connection is calculated (and stored for runtime)
# And the returned connection a's return to that.
function _fiber_connection_a(e::Edge, i::Int64, V::GKM_vector_bundle)
  return get_attribute(V, :connectionA)[(e, i)]
end

# Return the line bundle $\mathcal{O}_{\mathbb{P}^n}(d)$.
# TODO: document and export this, and maybe chose better way of linearization. 
function line_bundle_O(n::Int64, d::Int64)
  Pn = projective_space(GKM_graph, n)
  g = gens(Pn.M)
  GMtoM = ModuleHomomorphism(Pn.M, Pn.M, [g[i] for i in 1:n+1]);
  V = line_bundle(Pn, Pn.M, GMtoM, vcat([g[1]], [g[1] - d*_w(Pn, Edge(1, v)) for v in 2:n+1]))
  return V
end

# Return the line bundle $\mathcal{O}_{\mathbb{P}^n}(d)$.
# TODO: document and export this, and maybe chose better way of linearization. 
@doc raw"""
    vector_bundle_O(n::Int64, d::Vector{Int64})

Given $d=(d_1, d_2,\dots,d_r)$, return the vector bundle
$\mathcal{O}(d_1)\oplus\cdots\oplus\mathcal{O}(d_r)\rightarrow \mathbb{P}^n$,
linearized so that each summand has its own equivariant parameter.

# Example
Let us see $\mathcal{O}(-3)\oplus\mathcal{O}(0)\oplus\mathcal{O}(5)$ on $\mathbb{P}^3$.

```jldoctest vec_bdle_O_test
julia> V = GKMtools.vector_bundle_O(3, [-3, 0, 5])
GKM vector bundle of rank 3 over GKM graph with 4 nodes and valency 3 with weights:
1: (0, 0, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 1)
2: (3, -3, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1, 0), (-5, 5, 0, 0, 0, 0, 1)
3: (3, 0, -3, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1, 0), (-5, 0, 5, 0, 0, 0, 1)
4: (3, 0, 0, -3, 1, 0, 0), (0, 0, 0, 0, 0, 1, 0), (-5, 0, 0, 5, 0, 0, 1)
```

The acting torus has rank 7. The first 4 copies of $\mathbb{C}^\times$ act on
the base $\mathbb{P}^3$, the remaining 3 are just there to scale the fibres.
"""
function vector_bundle_O(n::Int64, d::Vector{Int64})
  r = length(d)
  Pn = enlarge_torus(projective_space(GKM_graph, n), r)
  g = gens(Pn.M)
  GMtoM = ModuleHomomorphism(Pn.M, Pn.M, [g[i] for i in 1:n+1+r]);
  # line_bdles = Vector{GKM_vector_bundle}()
  line_bdles = Vector{GKM_vector_bundle}(undef, r)
  ctr = 0
  for a in d
    ctr += 1
    L = line_bundle(Pn, Pn.M, GMtoM, vcat([g[ctr+n+1]], [g[ctr+n+1] - a*_w(Pn, Edge(1, v)) for v in 2:n+1]))
    # push!(line_bdles, L)
    line_bdles[ctr] = L
  end
  V = GKMtools.direct_sum(line_bdles...,)
  return V
end



@doc raw"""
    wedge_product(V::GKM_vector_bundle, n::Int64) -> GKM_vector_bundle

Return the wedge product, or external product, $\wedge^n V$.

# Example
Let us compute $\wedge^3 Q$ where `Q` is the universal quotient bundle of the Grassmannian $G(2, 5)$.
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKM_graph, 2, 5);

julia> wedge_product(Q, 2)
GKM vector bundle of rank 3 over GKM graph with 10 nodes and valency 6 with weights:
12: (0, 0, -1, -1, 0), (0, 0, -1, 0, -1), (0, 0, 0, -1, -1)
13: (0, -1, 0, -1, 0), (0, -1, 0, 0, -1), (0, 0, 0, -1, -1)
14: (0, -1, -1, 0, 0), (0, -1, 0, 0, -1), (0, 0, -1, 0, -1)
15: (0, -1, -1, 0, 0), (0, -1, 0, -1, 0), (0, 0, -1, -1, 0)
23: (-1, 0, 0, -1, 0), (-1, 0, 0, 0, -1), (0, 0, 0, -1, -1)
24: (-1, 0, -1, 0, 0), (-1, 0, 0, 0, -1), (0, 0, -1, 0, -1)
25: (-1, 0, -1, 0, 0), (-1, 0, 0, -1, 0), (0, 0, -1, -1, 0)
34: (-1, -1, 0, 0, 0), (-1, 0, 0, 0, -1), (0, -1, 0, 0, -1)
35: (-1, -1, 0, 0, 0), (-1, 0, 0, -1, 0), (0, -1, 0, -1, 0)
45: (-1, -1, 0, 0, 0), (-1, 0, -1, 0, 0), (0, -1, -1, 0, 0)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function wedge_product(V::GKM_vector_bundle, n::Int64)::GKM_vector_bundle
  @req n<=rank(V) "n is greater than the rank"

  return _wedge_and_sym_product(V, n, true)
end
@doc raw"""
    sym_product(V::GKM_vector_bundle, n::Int64) -> GKM_vector_bundle

Return the symmetric product $\mathrm{Sym}^n V$.

# Example
Let us compute $\mathrm{Sym}^2 Q$ where `Q` is the universal quotient bundle of the Grassmannian $G(3, 5)$.
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKM_graph, 3, 5);

julia> sym_product(Q, 2)
GKM vector bundle of rank 3 over GKM graph with 10 nodes and valency 6 with weights:
123: (0, 0, 0, -2, 0), (0, 0, 0, -1, -1), (0, 0, 0, 0, -2)
124: (0, 0, -2, 0, 0), (0, 0, -1, 0, -1), (0, 0, 0, 0, -2)
125: (0, 0, -2, 0, 0), (0, 0, -1, -1, 0), (0, 0, 0, -2, 0)
134: (0, -2, 0, 0, 0), (0, -1, 0, 0, -1), (0, 0, 0, 0, -2)
135: (0, -2, 0, 0, 0), (0, -1, 0, -1, 0), (0, 0, 0, -2, 0)
145: (0, -2, 0, 0, 0), (0, -1, -1, 0, 0), (0, 0, -2, 0, 0)
234: (-2, 0, 0, 0, 0), (-1, 0, 0, 0, -1), (0, 0, 0, 0, -2)
235: (-2, 0, 0, 0, 0), (-1, 0, 0, -1, 0), (0, 0, 0, -2, 0)
245: (-2, 0, 0, 0, 0), (-1, 0, -1, 0, 0), (0, 0, -2, 0, 0)
345: (-2, 0, 0, 0, 0), (-1, -1, 0, 0, 0), (0, -2, 0, 0, 0)
```
"""
function sym_product(V::GKM_vector_bundle, n::Int64)::GKM_vector_bundle

  return _wedge_and_sym_product(V, n, false)
end

function _wedge_and_sym_product(V::GKM_vector_bundle, n::Int64, wedged::Bool)::GKM_vector_bundle
  
  @req n>-1 "n must be non negative"
  
  if n == 1
    return V
  elseif n == 0
    return _zero_line_bundle(V)
  end

  G = V.gkm
  indices = wedged ? collect(powerset([i for i in 1:rank(V)], n, n)) : collect(with_replacement_combinations([i for i in 1:rank(V)], n))
  # rank_w = wedged ? binomial(rank(V), n) : binomial(rank(V) + n - 1, n)
  rank_w = length(indices)
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(_get_weight_type(G))}}(undef, n_vertices(G.g), rank_w)

  for v in 1:n_vertices(G.g)
    for r in 1:rank_w
      weightMatrix[v, r] = sum(i -> V.w[v, i], indices[r])
    end
  end

  res = vector_bundle(G, V.M, V.GMtoM, weightMatrix; calculateConnection = false)

  # If V has a connection, compute the induced connection on the wedge/symmetric product
  con_V = get_connection(V)
  if !isnothing(con_V)
    bundle_con = Dict{Tuple{Edge, Int64}, Int64}()

    for e in edges(G.g)
      for r in 1:rank_w
        # indices[r] is the r-th multi-index (set or multiset depending on wedged)
        # Apply the connection to each component
        dst_indices = sort([con_V[(e, i)] for i in indices[r]])

        # Find which index in indices corresponds to dst_indices
        dst_r = findfirst(idx -> sort(idx) == dst_indices, indices)
        @req !isnothing(dst_r) "Connection on wedge/symmetric product is not well-defined"

        bundle_con[(e, r)] = dst_r
        bundle_con[(reverse(e), dst_r)] = r
      end
    end

    res.con = bundle_con
    @req isvalid(bundle_con, res) "Wedge or sym product resulted in invalid bundle connection."
  end

  return res
end

@doc raw"""
    ^(V::GKM_vector_bundle, n::Number) -> GKM_vector_bundle

Return the tensor product $V^{\otimes n}$.

# Example
Let us compute the line bundle $l=\mathcal{O}(-4)$ of the Grassmannian $G(2, 4)$.
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKM_graph, 2, 4);

julia> S
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

  # @req rank(V) == 1 "currently tensor product implemented only for line bundles"

  if n == 1
    return V
  elseif n == 0
    return _zero_line_bundle(V)
  elseif rank(V) > 1
    return prod(i -> V, 1:n)
  end

  # For line bundles (rank 1), the n-th tensor power has a simple connection:
  # if V has connection con_V, then V^n has connection con_V (same connection)
  res = vector_bundle(V.gkm, V.M, V.GMtoM, n*V.w; calculateConnection = false)

  con_V = get_connection(V)
  if !isnothing(con_V)
    # For a line bundle, the connection on V^n is the same as the connection on V
    res.con = con_V
    @req isvalid(con_V, res) "Tensor power resulted in invalid bundle connection."
  end

  return res
end

@doc raw"""
    *(V::GKM_vector_bundle, W::GKM_vector_bundle) -> GKM_vector_bundle

Return the tensor product of `V` and `W`, that is $V \otimes W$.

# Example
Let us compute the vector bundle $S\otimes\mathcal{O}(-1)$ of the Grassmannian $G(2, 4)$.
```jldoctest
julia> S, Q = tautological_and_univ_bd(GKM_graph, 2, 4);

julia> S
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

julia> S * O_minus_one
GKM vector bundle of rank 2 over GKM graph with 6 nodes and valency 4 with weights:
12: (-2, -1, 0, 0), (-1, -2, 0, 0)
13: (-2, 0, -1, 0), (-1, 0, -2, 0)
14: (-2, 0, 0, -1), (-1, 0, 0, -2)
23: (0, -2, -1, 0), (0, -1, -2, 0)
24: (0, -2, 0, -1), (0, -1, 0, -2)
34: (0, 0, -2, -1), (0, 0, -1, -2)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function *(V::GKM_vector_bundle, W::GKM_vector_bundle)::GKM_vector_bundle

  # @req rank(V) == 1 "currently tensor product implemented only for line bundles"

  nv = n_vertices(V.gkm.g)
  rank_product = rank(V)*rank(W)
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(_get_weight_type(V.gkm))}}(undef, nv, rank_product)

  for v in 1:nv
    for i in 1:rank(V)
      for j in 1:rank(W)
        weightMatrix[v, (j-1)*rank(V) + i] = V.w[v, i] + W.w[v, j]
      end
    end
  end

  res = vector_bundle(V.gkm, V.M, V.GMtoM, weightMatrix; calculateConnection = false)

  # If both V and W have connections, compute the induced connection on the tensor product
  con_V = get_connection(V)
  con_W = get_connection(W)
  if !isnothing(con_V) && !isnothing(con_W)
    bundle_con = Dict{Tuple{Edge, Int64}, Int64}()

    for e in edges(V.gkm.g)
      for i in 1:rank(V)
        for j in 1:rank(W)
          # The (i,j)-th fiber at src(e) connects to (con_V[(e,i)], con_W[(e,j)])-th fiber at dst(e)
          i_dst = con_V[(e, i)]
          j_dst = con_W[(e, j)]

          # The index in the tensor product is (j-1)*rank(V) + i
          idx_src = (j-1)*rank(V) + i
          idx_dst = (j_dst-1)*rank(V) + i_dst

          bundle_con[(e, idx_src)] = idx_dst
          bundle_con[(reverse(e), idx_dst)] = idx_src
        end
      end
    end

    res.con = bundle_con
    @req isvalid(bundle_con, res) "Product resulted in invalid bundle connection"
  end

  return res
end

function _zero_line_bundle(V::GKM_vector_bundle)

  nv = n_vertices(V.gkm.g)
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{typeof(_get_weight_type(V.gkm))}}(undef, nv, 1)

  fill!(weightMatrix, 0*V.w[1, 1])

  res = vector_bundle(V.gkm, V.M, V.GMtoM, weightMatrix; calculateConnection = false)

  # The trivial bundle has a trivial connection: each fiber flag connects to itself
  # This is well-defined if the base GKM graph has a connection
  bundle_con = Dict{Tuple{Edge, Int64}, Int64}()
  for e in edges(V.gkm.g)
    bundle_con[(e, 1)] = 1
    bundle_con[(reverse(e), 1)] = 1
  end
  res.con = bundle_con

  return res
end

@doc raw"""
    gkm_line_bundle_of_toric(V::ToricLineBundle) -> GKM_vector_bundle

Return the GKM line bundle supported on the toric line bundle `V`.

# Example
Let us compute the line bundle of bidegree `[1,0]` on the Hirzebruch surface $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}\oplus \mathcal{O}_{\mathbb{P}^1}(4))$.
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4);

julia> V = toric_line_bundle(F4, [1,0]);

julia> gkm_line_bundle_of_toric(V)
GKM vector bundle of rank 1 over GKM graph with 4 nodes and valency 2 with weights:
1: (0, 0, -1, 1, 0)
2: (1, 0, -1, 0, 0)
3: (1, 0, -1, 0, 0)
4: (0, 0, -1, 1, 0)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function gkm_line_bundle_of_toric(V::ToricLineBundle)
  return gkm_vector_bundle_of_toric([V])
end

function gkm_vector_bundle_of_toric(V::ToricLineBundle)
  return gkm_vector_bundle_of_toric([V])
end

@doc raw"""
    gkm_vector_bundle_of_toric(E::Vector{ToricLineBundle}) -> GKM_vector_bundle

Return the GKM vector bundle supported on the direct sum of the toric line bundles in the vector `E`.

# Example
Let us compute the direct sum of bundles of bidegree `[1,0]` and `[1,1]` on the Hirzebruch surface $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}\oplus \mathcal{O}_{\mathbb{P}^1}(4))$.
```jldoctest
julia> F4 = hirzebruch_surface(NormalToricVariety, 4);

julia> V1 = toric_line_bundle(F4, [1,0]);

julia> V2 = toric_line_bundle(F4, [1,1]);

julia> E = [V1, V2];

julia> gkm_vector_bundle_of_toric(E)
GKM vector bundle of rank 2 over GKM graph with 4 nodes and valency 2 with weights:
1: (0, 0, -1, 0, 1, 0), (0, 0, 0, -1, -3, 1)
2: (1, 0, -1, 0, 0, 0), (-3, 0, 0, -1, 0, 1)
3: (1, 0, -1, 0, 0, 0), (1, 1, 0, -1, 0, 0)
4: (0, 0, -1, 0, 1, 0), (0, 1, 0, -1, 1, 0)
```
!!! warning
    All constructions involving vector bundles of the package are under develpment and will be expanded in the future.

"""
function gkm_vector_bundle_of_toric(E::Vector{ToricLineBundle})
  base = toric_variety(E[1])
  @req is_projective(base) "toric variety must be projective"
  @req is_smooth(base) "toric variety must be smooth, non-smooth not supported yet"

  v = total_space(E...)
  len = length(maximal_cones(v))
  g = Graph{Undirected}(len)
  M = free_module(ZZ, n_rays(v))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}()
  
  for sigma1 in 1:(len-1)
    for sigma2 in (sigma1+1):len
      count(x -> x in rays(maximal_cones(v)[sigma1]), rays(maximal_cones(v)[sigma2])) != (dim(v) - 1) && continue
      add_edge!(g, sigma1, sigma2)
      ray1 = findfirst(r -> !(r in rays(maximal_cones(v)[sigma2])), rays(maximal_cones(v)[sigma1]))
      W[Edge(sigma2, sigma1)] = _omega(v, sigma1, ray1, M)
    end
  end
  
  # gkm_graph(g, ["$i" for i in 1:len], M, W)
  add_dim = length(E) # additional dimension
  weightMatrix = Matrix{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}(undef, len, add_dim)

  scalars = gens(M)
  ud = QQFieldElem[]
  
  
  for j in 1:add_dim
    ray = ray_vector([i == (dim(v) - add_dim + j) for i in 1:dim(v)])
    for n_SIGMA1 in 1:len
      SIGMA1 = maximal_cones(v)[n_SIGMA1]
      # SIGMA1 = maximal_cones(v)[n_SIGMA1 == 1 ? len : (n_SIGMA1 == len ? 1 : n_SIGMA1)]

      for pol_ray in rays(polarize(SIGMA1))
        dot(ray, pol_ray) == 0 && continue
        ud = lcm(denominator.(pol_ray)) * pol_ray
        break
      end

      weightMatrix[n_SIGMA1, j] = -sum(k -> Int64(dot(rays(v)[k], ud))*scalars[k], 1:n_rays(v))

    end
  end

  # return gkm_graph(g, ["$i" for i in 1:len], M, W)

  GMtoM = ModuleHomomorphism(M, M, [gens(M)[i] for i in 1:n_rays(v)])

  vector_bundle(gkm_graph(g, ["$i" for i in 1:len], M, W), M, GMtoM, weightMatrix; calculateConnection = true)
end