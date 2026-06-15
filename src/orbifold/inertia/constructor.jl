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
#                (b) the flag embedding at the (v->w) flag sends g to some
#                    element f in the flag isotropy group, AND
#                    the flag embedding at the (w->v) flag sends h to the
#                    same element f.
#              In other words g and h are "conjugate via the flag":
#              they both restrict to the same element of the flag isotropy.
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
    apply_embedding(embedding::ZZMatrix, g::Vector{Int},
                    src_structure::Vector{Int},
                    tgt_structure::Vector{Int}) -> Vector{Int}

Apply the linear map `embedding` (a k×m integer matrix) to the group element
`g ∈ Z/n1×...×Z/nm`, reducing each coordinate modulo the target group order.

`embedding` rows index the target group factors, columns the source factors.
"""
function apply_embedding(
  embedding::ZZMatrix,
  g::Vector{Int},
  src_structure::Vector{Int},
  tgt_structure::Vector{Int},
)::Vector{Int}
  # Nothing to do for trivial flag isotropy group
  isempty(tgt_structure) && return Int[]

  k = length(tgt_structure)   # target (flag isotropy) rank
  m = length(src_structure)   # source (vertex isotropy) rank
  @assert size(embedding, 1) == k && size(embedding, 2) == m """
     Embedding matrix size $(size(embedding)) is incompatible with
     source rank $m → target rank $k.
 """

  result = zeros(Int, k)
  for i in 1:k
    val = 0
    for j in 1:m
      val += Int(embedding[i, j]) * g[j]
    end
    result[i] = mod(val, tgt_structure[i])
  end
  return result
end

###############################################################################
# 2.  Constructing the inertia stack vertex / edge data
###############################################################################

"""
    InertiaVertex{V}

A vertex of the inertia stack: original vertex label plus a group element.
"""
struct InertiaVertex{V}
  original_vertex::V            # label in the original GKM graph
  vertex_index::Int             # index in the original graph (1-based)
  sector_element::Vector{Int}   # element of G_v (0 = untwisted sector)
end

Base.show(io::IO, iv::InertiaVertex) =
  print(io, "(", iv.original_vertex, ", ", iv.sector_element, ")")

"""
    inertia_stack(X::OrbifoldGKMGraph) -> OrbifoldGKMGraph

Construct the inertia stack IX of the orbifold GKM graph X.

The returned graph is again an OrbifoldGKMGraph whose:
  - vertices are InertiaVertex objects (pairs (v, g)),
  - edges / flags are inherited from matched flag pairs,
  - isotropy and flag-isotropy data come from the original vertex / flag data.
"""
function inertia_stack(X::OrbifoldGKMGraph{R,V,F}) where {R,V,F}
  core = X.core
  n = nv(core.g)               # number of original vertices
  orig_labels = core.labels           # V-typed labels

  ###########################################################################
  # Step 1: Enumerate all inertia vertices (v, g)
  ###########################################################################

  # For each original vertex v, list every element g ∈ G_v.
  # The isotropy group structure is X.vertex_isotropy[v].isotropy_group.

  inertia_vertices = InertiaVertex{V}[]   # flat list of all (v, g) pairs
  # For fast lookup: original_vertex_index, element → inertia_vertex_index
  vertex_to_inertia = Dict{Tuple{Int,Vector{Int}},Int}()

  for v in 1:n
    gv_structure = X.vertex_isotropy[v].isotropy_group  # e.g. [3, 2]
    for g in abelian_group_elements(gv_structure)
      iv = InertiaVertex{V}(orig_labels[v], v, g)
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
    fi_src_idx, fi_tgt_idx = core.edge_flags[e]
    # flag v→w: embedding from G_v into G_{flag}
    flag_vw = X.flag_isotropy[v][fi_src_idx]
    # flag w→v: embedding from G_w into G_{flag}
    flag_wv = X.flag_isotropy[w][fi_tgt_idx]

    gv_structure = X.vertex_isotropy[v].isotropy_group
    gw_structure = X.vertex_isotropy[w].isotropy_group
    gf_structure = flag_vw.isotropy_group   # flag isotropy group structure

    # For each g ∈ G_v, compute its image f_g in G_{flag}.
    # For each h ∈ G_w, compute its image f_h in G_{flag}.
    # An inertia edge (v,g)-(w,h) exists iff f_g == f_h.

    # Build a map: flag_element → list of g ∈ G_v mapping there
    fiber_v = Dict{Vector{Int},Vector{Vector{Int}}}()
    for g in abelian_group_elements(gv_structure)
      f = apply_embedding(flag_vw.embedding, g, gv_structure, gf_structure)
      push!(get!(fiber_v, f, Vector{Int}[]), g)
    end

    for h in abelian_group_elements(gw_structure)
      f = apply_embedding(flag_wv.embedding, h, gw_structure, gf_structure)
      if haskey(fiber_v, f)
        for g in fiber_v[f]
          i = vertex_to_inertia[(v, g)]
          j = vertex_to_inertia[(w, h)]
          # Store with i < j to avoid duplicate undirected edges
          push!(inertia_edges, i < j ? (i, j) : (j, i))
        end
      end
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

  # Inertia vertex labels are InertiaVertex objects.
  # If we want to stay in the same type universe we can use strings;
  # here we keep them as-is and parameterise the new graph on InertiaVertex{V}.
  inertia_labels = inertia_vertices   # Vector{InertiaVertex{V}}

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

    # Retrieve the original edge (order may vary)
    orig_edge = has_edge(core.g, v, w) ? Edge(v, w) : Edge(w, v)
    orig_fi_src, orig_fi_tgt = core.edge_flags[orig_edge]

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

    orig_edge = has_edge(core.g, v, w) ? Edge(v, w) : Edge(w, v)
    orig_fi_src, orig_fi_tgt = core.edge_flags[orig_edge]

    # The inertia flag isotropy is the same as the original flag isotropy.
    push!(inertia_flag_isotropy[i], X.flag_isotropy[v][orig_fi_src])
    push!(inertia_flag_isotropy[j], X.flag_isotropy[w][orig_fi_tgt])
  end

  ###########################################################################
  # Step 6: Assemble the OrbifoldGKMGraph
  ###########################################################################

  inertia_core = GKMCombinatorialData{R,InertiaVertex{V},F}(
    g_inertia,
    core.M,
    inertia_labels,
    inertia_flags,
    inertia_edge_flags,
  )

  return OrbifoldGKMGraph{R,InertiaVertex{V},F}(
    inertia_core,
    inertia_vertex_isotropy,
    inertia_flag_isotropy,
  )
end

###############################################################################
# 3.  Convenience: list the twisted sectors
###############################################################################

"""
    twisted_sectors(IX::OrbifoldGKMGraph{R, InertiaVertex{V}, F})
        -> Vector{InertiaVertex{V}}

Return all inertia vertices belonging to a twisted sector,
i.e. those whose sector element is not the identity (zero vector).
"""
function twisted_sectors(IX::OrbifoldGKMGraph{R,InertiaVertex{V},F}) where {R,V,F}
  return filter(iv -> any(!iszero, iv.sector_element), IX.core.labels)
end

"""
    untwisted_sector(IX::OrbifoldGKMGraph{R, InertiaVertex{V}, F})
        -> Vector{InertiaVertex{V}}

Return all inertia vertices belonging to the untwisted sector
(sector element = zero vector).
"""
function untwisted_sector(IX::OrbifoldGKMGraph{R,InertiaVertex{V},F}) where {R,V,F}
  return filter(iv -> all(iszero, iv.sector_element), IX.core.labels)
end

"""
    sector_count(IX::OrbifoldGKMGraph{R, InertiaVertex{V}, F}) -> Int

Number of distinct twisted sectors (not counting the untwisted one).
"""
function sector_count(IX::OrbifoldGKMGraph{R,InertiaVertex{V},F}) where {R,V,F}
  return length(twisted_sectors(IX))
end
