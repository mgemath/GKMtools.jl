@doc raw"""
    serialize_gkm_graph(path, G)

Store in a file the GKM graph `G` in the path given by the string `path`.
The cohomology of the graph is not stored, but it is recomputed when loading the graph. Use `deserialize_gkm_graph(path)` to restore the graph, with Oscar and GKMtools loaded.

!!! warning "Warning"
    Currently, only generalized flag graphs are supported. 

!!! warning "Warning"
    The file format could be changed in future versions of GKMtools, so retro-compatibility is not guaranteed.

```jldoctest
julia> G = generalized_gkm_flag(root_system(:A, 3), [1, 2]);

julia> path = tempname();

julia> serialize_gkm_graph(path, G);

julia> G = deserialize_gkm_graph(path)
GKM graph with 4 nodes, valency 3 and axial function:
s3 -> id => (0, 0, 1, -1)
s2*s3 -> id => (0, 1, 0, -1)
s2*s3 -> s3 => (0, 1, -1, 0)
s1*s2*s3 -> id => (1, 0, 0, -1)
s1*s2*s3 -> s3 => (1, 0, -1, 0)
s1*s2*s3 -> s2*s3 => (1, -1, 0, 0)
Birkhoff-Grothendieck connection for GKM graph with 4 nodes and valency 3

julia> rm(path);
```
"""
function serialize_gkm_graph(path::AbstractString,
    G::GKMGraph{C,GeneralizedFlagVertex,FlagWeight{C}}) where C
  R = root_system(parent(flag(first(G.core.labels))))
  native(x) = Rational{BigInt}(x)
  edgekey(e) = (src(e), dst(e))
  h = G.H2
  data = (
    format = "GKMtools.generalized_flag", version = 1,
    cartan = BigInt.(Matrix(cartan_matrix(R))),
    rational = C == QQFieldElem, lattice_rank = rank(G.core.M),
    labels = get_string.(G.core.labels),
    words = [Int.(word(flag(v))) for v in G.core.labels],
    weights = [[native.([f.weight[i] for i in 1:rank(G.core.M)]) for f in fs] for fs in G.core.flags],
    edge_flags = Dict(edgekey(e) => ij for (e, ij) in G.core.edge_flags),
    transport = Dict(edgekey(e) => copy(p) for (e, p) in G.connection.transport),
    coefficients = Dict(edgekey(e) => native.(a) for (e, a) in G.connection.a),
    connection_type = G.connection.connection_type,
    h2_rank = rank(h.H2),
    edge_to_gen = Dict(edgekey(e) => i for (e, i) in h.edge_to_gen),
    quotient = [BigInt.([h.quotient(g)[i] for i in 1:rank(h.H2)]) for g in gens(h.edge_lattice)],
    chern = [BigInt(h.chern(g)[1]) for g in gens(h.H2)],
    ray_sum = native.(h.ray_sum),
  )
  temporary, io = mktemp(dirname(abspath(path)))
  try
    Serialization.serialize(io, data)
    close(io)
    mv(temporary, path; force=true)
  finally
    isopen(io) && close(io)
    isfile(temporary) && rm(temporary)
  end
  return G
end

@doc raw"""
    deserialize_gkm_graph(path)

Restore a generalized GKM graph saved by [`serialize_gkm_graph`](@ref).
"""
function deserialize_gkm_graph(path::AbstractString)
  d = Serialization.deserialize(path)
  (d isa NamedTuple && get(d, :format, nothing) == "GKMtools.generalized_flag" &&
    get(d, :version, nothing) == 1) || throw(ArgumentError("Unsupported GKM graph file format or version"))
  R = root_system(matrix(ZZ, d.cartan))
  W = weyl_group(R)
  C = d.rational ? QQFieldElem : ZZRingElem
  ring = d.rational ? QQ : ZZ
  M = free_module(ring, d.lattice_rank)
  labels = [GeneralizedFlagVertex(l, W(w)) for (l, w) in zip(d.labels, d.words)]
  fs = [FlagWeight{C}[FlagWeight{C}(M(ring.(w))) for w in ws] for ws in d.weights]
  g = Graph{Undirected}(length(labels))
  ef = Dict{Edge,Tuple{Int,Int}}()
  for ((s, t), ij) in d.edge_flags
    add_edge!(g, s, t)
    ef[Edge(s, t)] = ij
  end
  core = GKMCombinatorialData{C,GeneralizedFlagVertex,FlagWeight{C}}(g, M, labels, fs, ef)
  con = Connection{C}(
    Dict(Edge(s, t) => p for ((s, t), p) in d.transport),
    Dict(Edge(s, t) => ring.(a) for ((s, t), a) in d.coefficients),
    d.connection_type,
  )
  h = if isempty(d.edge_to_gen)
    _zero_GKM_H2()
  else
    E = free_module(ZZ, length(d.quotient))
    H = free_module(ZZ, d.h2_rank)
    Z = free_module(ZZ, 1)
    quotient = ModuleHomomorphism(E, H, [H(ZZ.(v)) for v in d.quotient])
    chern = ModuleHomomorphism(H, Z, [Z([ZZ(c)]) for c in d.chern])
    GKM_H2(E, H, Dict(Edge(s, t) => i for ((s, t), i) in d.edge_to_gen),
      quotient, cone_from_inequalities(-identity_matrix(QQ, d.h2_rank)), QQ.(d.ray_sum), chern)
  end
  return GKMGraph{C,GeneralizedFlagVertex,FlagWeight{C}}(
    core, con, create_cohomology(rank(M), length(labels)), h, nothing)
end
