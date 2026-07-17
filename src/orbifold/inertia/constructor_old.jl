###############################################################################
# Inertia Stack of an Orbifold GKM Graph
#
# The inertia stack IXof an orbifold GKM graph X has:
#
#   Vertices : pairs (v, g) where v is a vertex of X and g is an element
#              of the isotropy group G_v at v.  The identity element g = 0
#              gives the untwisted sector; all other g give twisted sectors.
#
#   Edges    : (v, g) -- (w, h)  iff
#                (a) there is an edge v--w in X, AND
#                (b) there is an element f in the flag isotropy group whose
#                    images under H_flag -> G_v and H_flag -> G_w are g and h.
#              In other words, the sector element must extend over the
#              one-dimensional stratum represented by the flag.
#
#   Isotropy at (v, g) : the full G_v (abelian, so centralizer = G_v).
#
#   Flag isotropy at (v,g)->(w,h) : the flag isotropy at the original flag v->w.
#
#   Flag weights : inherited unchanged from the original edge.
###############################################################################

###############################################################################
# 1.  Helpers for finite abelian group arithmetic
###############################################################################

"""
    abelian_group_elements(structure::Vector{Int}) -> Vector{Vector{Int}}

Enumerate all elements of Z/n1 × Z/n2 × ... × Z/nk.
Returns a vector of tuples represented as `Vector{Int}`.
The identity is the zero vector.
"""
function abelian_group_elements(structure::Vector{Int})::Vector{Vector{Int}}
  isempty(structure) && return [Int[]]          # trivial group: one element
  result = [Int[]]
  for n in structure
    result = [vcat(prefix, [i]) for prefix in result for i in 0:(n - 1)]
  end
  return result
end

"""
    apply_embedding(embedding::ZZMatrix, h::Vector{Int},
                    src_structure::Vector{Int},
                    tgt_structure::Vector{Int}) -> Vector{Int}

Apply the homomorphism encoded by `embedding` to a flag isotropy element.
Here `embedding` is the map H_flag -> G_vertex: rows index the target
vertex isotropy factors and columns index the source flag isotropy factors.
"""
function apply_embedding(
  embedding::ZZMatrix,
  h::Vector{Int},
  src_structure::Vector{Int},
  tgt_structure::Vector{Int},
)::Vector{Int}
  k = length(tgt_structure)
  m = length(src_structure)
  m == 0 && return zeros(Int, k)

  @assert size(embedding, 1) == k && size(embedding, 2) == m """
     Embedding matrix size $(size(embedding)) is incompatible with
     source rank $m -> target rank $k.
 """

  result = zeros(Int, k)
  for i in 1:k
    val = 0
    for j in 1:m
      val += Int(embedding[i, j]) * h[j]
    end
    result[i] = mod(val, tgt_structure[i])
  end
  return result
end

function sector_age(I::OrbifoldVertexIsotropy, g::Vector{Int})
  age = QQ(0)
  for j in 1:size(I.tangent_rep, 2)
    exponent = QQ(0)
    for i in 1:length(I.isotropy_group)
      exponent += QQ(g[i] * Int(I.tangent_rep[i, j]), I.isotropy_group[i])
    end
    age += exponent - floor(exponent)
  end
  return age
end

###############################################################################
# 2.  Constructing the inertia stack vertex / edge data
###############################################################################

"""
    InertiaStackVertex{V}

A vertex of the inertia stack: original vertex label plus a group element.
"""
struct InertiaStackVertex{V} <: AbstractVertex
  label::String
  original_vertex::V            # label in the original GKM graph
  vertex_index::Int             # index in the original graph (1-based)
  sector_element::Vector{Int}   # element of G_v (0 = untwisted sector)
  age::QQFieldElem
end

function InertiaStackVertex(
  original_vertex::V,
  vertex_index::Int,
  sector_element::Vector{Int},
  age::QQFieldElem = QQ(0),
) where {V}
  label = "($(_vertex_label(original_vertex)), $(sector_element))"
  return InertiaStackVertex{V}(label, original_vertex, vertex_index, sector_element, age)
end

_vertex_label(v::AbstractVertex) = get_string(v)
_vertex_label(v) = string(v)

Base.show(io::IO, iv::InertiaStackVertex) =
  print(io, iv.label)

function _edge_flag_indices(core::GKMCombinatorialData, v::Int, w::Int)
  e = Edge(v, w)
  if haskey(core.edge_flags, e)
    return core.edge_flags[e]
  end

  reverse_e = Edge(w, v)
  if haskey(core.edge_flags, reverse_e)
    reverse_src, reverse_dst = core.edge_flags[reverse_e]
    return reverse_dst, reverse_src
  end

  error("Edge $v-$w not found in edge_flags")
end

"""
    InertiaStack{R,V,F,X}

The inertia stack of an orbifold GKM graph.

It stores the original orbifold graph in `base` and the GKM graph of the
inertia stack in the usual `core`, `vertex_isotropy`, and `flag_isotropy`
fields, so the standard `AbstractGKMGraph` interface applies directly.
"""
struct InertiaStack{
  R,
  V,
  F,
  X<:AbstractOrbifoldGKMGraph{R,V,F},
} <: AbstractOrbifoldGKMGraph{R,InertiaStackVertex{V},F}
  base::X
  core::GKMCombinatorialData{R,InertiaStackVertex{V},F}
  vertex_isotropy::Vector{OrbifoldVertexIsotropy}
  flag_isotropy::Vector{Vector{OrbifoldFlagIsotropy}}
end

function Base.show(io::IO, IX::InertiaStack)
  print(io, "InertiaStack with $(num_vertices(IX)) stacky vertices")
end

function Base.show(io::IO, ::MIME"text/plain", IX::InertiaStack)
  print(
    io,
    "Inertia stack of an orbifold GKM graph with $(num_vertices(IX)) nodes, $(num_edges(IX)) edges, and $(sector_count(IX)) twisted vertices",
  )
  for e in edges(IX)
    print(io, "\n$(label(IX, src(e))) -> $(label(IX, dst(e))) => $(weight(IX, e))")
  end
end

"""
    inertia_stack(X::OrbifoldGKMGraph) -> InertiaStack

Construct the inertia stack IX of the orbifold GKM graph X.

The returned object is an `InertiaStack` whose:
  - vertices are InertiaStackVertex objects (pairs (v, g)),
  - edges / flags are inherited from matched flag pairs,
  - isotropy and flag-isotropy data come from the original vertex / flag data.
"""
function inertia_stack(X::AbstractOrbifoldGKMGraph{R,V,F}) where {R,V,F}
  core = X.core
  n = nv(core.g)               # number of original vertices
  orig_labels = core.labels           # V-typed labels

  ###########################################################################
  # Step 1: Enumerate all inertia vertices (v, g)
  ###########################################################################

  # For each original vertex v, list every element g ∈ G_v.
  # The isotropy group structure is X.vertex_isotropy[v].isotropy_group.

  inertia_vertices = InertiaStackVertex{V}[]   # flat list of all (v, g) pairs
  # For fast lookup: original_vertex_index, element → inertia_vertex_index
  vertex_to_inertia = Dict{Tuple{Int,Vector{Int}},Int}()

  for v in 1:n
    gv_structure = X.vertex_isotropy[v].isotropy_group  # e.g. [3, 2]
    for g in abelian_group_elements(gv_structure)
      iv = InertiaStackVertex(orig_labels[v], v, g, sector_age(X.vertex_isotropy[v], g))
      push!(inertia_vertices, iv)
      vertex_to_inertia[(v, g)] = length(inertia_vertices)
    end
  end

  N = length(inertia_vertices)   # total number of inertia vertices

  ###########################################################################
  # Step 2: Build the inertia graph — determine which (v,g)-(w,h) pairs
  #         are connected.
  ###########################################################################

  # We will collect edges as (inertia_src_idx, inertia_tgt_idx).
  # Because Oscar's Graph is undirected we add each pair once.

  inertia_edges = Tuple{Int,Int}[]   # (i, j) with i < j

  for e in edges(core.g)
    v = src(e)
    w = dst(e)

    # Retrieve flag isotropy data for both directed flags on this edge.
    fi_src_idx, fi_tgt_idx = _edge_flag_indices(core, v, w)
    # flag v→w: embedding from G_v into G_{flag}
    flag_vw = X.flag_isotropy[v][fi_src_idx]
    # flag w→v: embedding from G_w into G_{flag}
    flag_wv = X.flag_isotropy[w][fi_tgt_idx]

    gv_structure = X.vertex_isotropy[v].isotropy_group
    gw_structure = X.vertex_isotropy[w].isotropy_group
    gf_structure = flag_vw.isotropy_group   # flag isotropy group structure

    # The flag embedding is H_flag -> G_vertex. A sector edge is
    # obtained by choosing a flag isotropy element and mapping it to both
    # endpoint stabilizers. If the flag isotropy is trivial, this connects
    # only the identity sectors.
    for f in abelian_group_elements(gf_structure)
      g = apply_embedding(flag_vw.embedding, f, gf_structure, gv_structure)
      h = apply_embedding(flag_wv.embedding, f, gf_structure, gw_structure)
      i = vertex_to_inertia[(v, g)]
      j = vertex_to_inertia[(w, h)]
      push!(inertia_edges, i < j ? (i, j) : (j, i))
    end
  end

  # Deduplicate (shouldn't happen in practice, but be safe)
  unique!(inertia_edges)

  ###########################################################################
  # Step 3: Build Oscar graph and flag data
  ###########################################################################

  g_inertia = Graph{Undirected}(N)
  for (i, j) in inertia_edges
    add_edge!(g_inertia, i, j)
  end

  # Inertia vertex labels are InertiaStackVertex objects.
  # If we want to stay in the same type universe we can use strings;
  # here we keep them as-is and parameterise the new graph on InertiaStackVertex{V}.
  inertia_labels = inertia_vertices   # Vector{InertiaStackVertex{V}}

  ###########################################################################
  # Step 4: Build flags for the inertia graph
  #
  # For each inertia edge (v,g)-(w,h), the flag weight is inherited from
  # the original edge v-w (same F-typed weight).
  ###########################################################################

  # We need to build:
  #   flags::Vector{Vector{F}}          — for each inertia vertex
  #   edge_flags::Dict{Edge,Tuple{Int,Int}}
  #
  # Strategy: for each inertia edge in order, look up the original edge,
  # read the flag weights, and append to the flag lists of the two endpoints.

  inertia_flags = [Vector{F}() for _ in 1:N]
  inertia_edge_flags = Dict{Edge,Tuple{Int,Int}}()

  for (i, j) in inertia_edges
    iv = inertia_vertices[i]
    jv = inertia_vertices[j]
    v = iv.vertex_index
    w = jv.vertex_index

    orig_fi_src, orig_fi_tgt = _edge_flag_indices(core, v, w)

    # Flag weight at source v (going towards w)
    fw_src = core.flags[v][orig_fi_src]
    # Flag weight at target w (going towards v)
    fw_tgt = core.flags[w][orig_fi_tgt]

    e_inertia = Edge(i, j)

    push!(inertia_flags[i], fw_src)
    src_fi = length(inertia_flags[i])

    push!(inertia_flags[j], fw_tgt)
    tgt_fi = length(inertia_flags[j])

    inertia_edge_flags[e_inertia] = (src_fi, tgt_fi)
  end

  ###########################################################################
  # Step 5: Build isotropy data for inertia vertices and flags
  #
  # Isotropy at inertia vertex (v, g) : full G_v (abelian ⟹ centralizer = G_v).
  # Isotropy at inertia flag (v,g)→(w,h): flag isotropy at the original flag.
  ###########################################################################

  inertia_vertex_isotropy = OrbifoldVertexIsotropy[]
  for iv in inertia_vertices
    # Centraliser of g in abelian G_v is G_v itself.
    push!(inertia_vertex_isotropy, X.vertex_isotropy[iv.vertex_index])
  end

  inertia_flag_isotropy = [Vector{OrbifoldFlagIsotropy}() for _ in 1:N]
  for (i, j) in inertia_edges
    iv = inertia_vertices[i]
    jv = inertia_vertices[j]
    v = iv.vertex_index
    w = jv.vertex_index

    orig_fi_src, orig_fi_tgt = _edge_flag_indices(core, v, w)

    # The inertia flag isotropy is the same as the original flag isotropy.
    push!(inertia_flag_isotropy[i], X.flag_isotropy[v][orig_fi_src])
    push!(inertia_flag_isotropy[j], X.flag_isotropy[w][orig_fi_tgt])
  end

  ###########################################################################
  # Step 6: Assemble the OrbifoldGKMGraph
  ###########################################################################

  inertia_core = GKMCombinatorialData{R,InertiaStackVertex{V},F}(
    g_inertia,
    core.M,
    inertia_labels,
    inertia_flags,
    inertia_edge_flags,
  )

  return InertiaStack{R,V,F,typeof(X)}(
    X,
    inertia_core,
    inertia_vertex_isotropy,
    inertia_flag_isotropy,
  )
end

###############################################################################
# 3.  Convenience: list the twisted sectors
###############################################################################

"""
    twisted_sectors(IX::InertiaStack{R,V,F})
        -> Vector{InertiaStackVertex{V}}

Return all inertia vertices belonging to a twisted sector,
i.e. those whose sector element is not the identity (zero vector).
"""
function twisted_sectors(IX::InertiaStack{R,V,F}) where {R,V,F}
  return filter(iv -> any(!iszero, iv.sector_element), IX.core.labels)
end

"""
    untwisted_sector(IX::InertiaStack{R,V,F})
        -> Vector{InertiaStackVertex{V}}

Return all inertia vertices belonging to the untwisted sector
(sector element = zero vector).
"""
function untwisted_sector(IX::InertiaStack{R,V,F}) where {R,V,F}
  return filter(iv -> all(iszero, iv.sector_element), IX.core.labels)
end

"""
    sector_count(IX::InertiaStack{R,V,F}) -> Int

Number of twisted-sector vertices (not counting the untwisted sector vertices).
"""
function sector_count(IX::InertiaStack{R,V,F}) where {R,V,F}
  return length(twisted_sectors(IX))
end
