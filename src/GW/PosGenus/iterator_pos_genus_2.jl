#########Functions relative to the graphs#########

function compute_threshold_for_progress_bar(max_genus::Int64, max_edges::Int64, n_marks, nc)::Int64

  # genus 0 case
  number_trees = A000055(max_edges + 1)
  threshold = (length(keys(nc))) * sum(vert -> number_trees[vert] * sum(genus -> length(weak_compositions(genus, vert)), 0:max_genus) * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:(max_edges + 1))

  max_genus == 0 && return threshold

  for (top_genus, n_vert) in Iterators.product(1:max_genus, 2:(max_edges + 1)) # we fix the top_genus and the number of vertices

    n_vert + top_genus - 1 > max_edges && continue # respect max_edges
    (top_genus > 0) && n_vert < ceil(Int64, (3 + sqrt(1 + 8*top_genus)) / 2) && continue # skip impossible cases
    n_den_dist_and_marks = sum(genus -> length(weak_compositions(genus, n_vert)), 0:(max_genus-top_genus)) * binomial(n_vert+n_marks-1, n_marks)

    for g6 in my_geng(top_genus, n_vert) # generation of graphs

      M = graph6_to_adjacency_matrix(g6)
      top_graph_Graphs = Graphs.SimpleGraph(M) # graph in Graphs of the current iteration
      threshold += count_homomorphisms(top_graph_Graphs, nc) * n_den_dist_and_marks
      
    end
  end

  return threshold
end

# function count_homomorphisms(G::SimpleGraph, H::SimpleGraph)
function count_homomorphisms(G::SimpleGraph, adjH)
  nG = Graphs.nv(G)
  # nH = nv(H)
  nH = length(keys(adjH))
  adjG = Graphs.SimpleGraphs.adj(G)
  # adjH = Graphs.SimpleGraphs.adj(H)

  # order vertices of G by descending degree
  order = sort(1:nG, by = v -> -Graphs.degree(G, v))
  pos = zeros(Int, nG)
  for (i, v) in enumerate(order)
    pos[v] = i
  end

  assignment = zeros(Int, nG)

  function backtrack(i)
    i > nG && return 1

    v = order[i]
    total = 0

    for h in 1:nH
      ok = true
      for u in adjG[v]
        if pos[u] < i && !(assignment[u] in adjH[h])
          ok = false
          break
        end
      end
      if ok
        assignment[v] = h
        total += backtrack(i + 1)
      end
    end

    return total
  end

  return backtrack(1)
end

function my_geng(g::Int64, n::Int64)::Base.EachLine{IOBuffer}

  # (g > 0) && n < ceil(Int, (3 + sqrt(1 + 8*g)) / 2) && return eachline(IOBuffer("")) # skip impossible cases, already checked in main loop
  n_edge = g + n - 1 # number of edges for connected graph with genus g and n vertices
  cmd = nauty_jll.geng_path * " -c -q $n $n_edge:$n_edge"
  iter_of_g6 = eachline(IOBuffer(read(`sh -c $cmd`)))

  return iter_of_g6
end

function compute_aut(top_graph_Graphs::Graphs.SimpleGraph{Int64})::Int64
  return Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs)
end

function compute_internal_aut_multiedges(multiedges::Tuple{Vararg{Vector{Int64}}})::Int64

  ans = 1
  for multi_edge in multiedges

    for m in unique(multi_edge)
      ans *= factorial(count(==(m), multi_edge))
    end

  end
  return ans
end


function graph6_to_adjacency_matrix(s::String)::Matrix{Bool}
  bytes = Vector{UInt8}(s)
  n = 0
  offset = 0

  if bytes[1] != 126
    n = bytes[1] - 63
    offset = 1
  else
    n = (bytes[2] - 63) << 12 | (bytes[3] - 63) << 6 | (bytes[4] - 63)
    offset = 4
  end

  total_edges = n * (n - 1) ÷ 2
  adj = zeros(Bool, n, n)

  if total_edges == 0
    return adj
  end

  bit_index = 0
  for i in offset+1:length(bytes)
    byte_val = bytes[i] - 63
    for j in 5:-1:0
      if bit_index < total_edges
        bit = (byte_val >> j) & 1
        if bit == 1
          # Calculate row and column from the bit index
          # The k-th bit corresponds to edge (i, j) where:
          # i = row, j = column, with i < j
          # The indexing is column-major for the upper triangle
          # Formula: j = floor((sqrt(8*k + 1) + 1)/2)
          #          i = k - j*(j-1)//2
          k = bit_index
          j_val = floor(Int, (sqrt(8*k + 1) + 1) / 2)
          i_val = k - j_val*(j_val - 1) ÷ 2
          adj[i_val+1, j_val+1] = Bool(1)
          adj[j_val+1, i_val+1] = Bool(1)
        end
        bit_index += 1
      end
    end
  end

  return adj
end



######## Generations of colorings
function collect_all_uniques_cols_W_C_N_flatmap(top_graph_Graphs::Graphs.SimpleGraph{Int64}, nc::Dict{Int64,Vector{Int64}}, top_aut::Int64)

  return Iterators.flatmap(comb -> unique_col_fixed_combination_w_counting_and_numering(top_graph_Graphs, nc, comb, top_aut), Combinatorics.with_replacement_combinations(1:length(nc), Graphs.nv(top_graph_Graphs)))   
end

function unique_col_fixed_combination_w_counting_and_numering(top_graph_Graphs::Graphs.SimpleGraph{Int64}, nc::Dict{Int64,Vector{Int64}}, comb::Vector{Int64}, top_aut::Int64)
  ans = Set{Tuple{Vector{Int64}, Int64}}()
  seen = Dict{Vector{Int64}, Int64}()

  total_number = length(Combinatorics.multiset_permutations(comb, length(comb)))

  for c in Combinatorics.multiset_permutations(comb, length(comb))

    all(e -> c[Graphs.src(e)] in nc[c[Graphs.dst(e)]], Graphs.edges(top_graph_Graphs)) || continue # check if coloring is valid
    
    found = false
    
    for color2 in ans
      seen[color2[1]] == 0 && continue

      color_rel(u, v) = (c[u] == color2[1][v])
      if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_rel)
        found = true
        seen[color2[1]] -= 1
        break
      end
    end

    if !found
      color_rel_2(u, v) = (c[u] == c[v])
      aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_rel_2)
      push!(ans, (c, aut))

      divis = div(top_aut, aut)
      seen[c] = divis - 1
      total_number -= divis
      total_number == 0 && break
    end
  end
  
  return ans
end

######## Generations of genus distributions

function genus_distribution_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, col, col_aut::Int64, max_genus::Int64, genus::Int64)

  ((col_aut == 1) || (max_genus == genus)) && return Iterators.zip(map(x -> collect(x), weak_compositions(max_genus - genus, Graphs.nv(top_graph_Graphs))), Iterators.cycle([col_aut]))

  return Iterators.flatmap(p -> unique_gen_dist_fixed_partition(top_graph_Graphs, p, col), Combinatorics.partitions(max_genus - genus))
end

function unique_gen_dist_fixed_partition(top_graph_Graphs::Graphs.SimpleGraph{Int64}, p::Vector{Int64}, col)
  
  ans = Set{Tuple{Vector{Int64}, Int64}}()
  length(p) > Graphs.nv(top_graph_Graphs) && return ans # impossible case

  extended_p = vcat(p, zeros(Int64, Graphs.nv(top_graph_Graphs) - length(p)))

  for perm_extended_p in Combinatorics.multiset_permutations(extended_p, Graphs.nv(top_graph_Graphs))

    found = false
    
    for (gen_dist, gen_dist_aut) in ans
      
      color_1(v, u) = (perm_extended_p[v] == gen_dist[u]) && (col[v] == col[u])
      
      if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_1) # check isomorphism
        found = true
        # seen[color2] -= 1 # use one copy of this coloring
        break
      end
    end

    if !found # new coloring if not found
      color_rel_2(u, v) = (perm_extended_p[v] == perm_extended_p[u]) && (col[v] == col[u])
      aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_rel_2) # count automorphisms of the coloring
      push!(ans, (perm_extended_p, aut))
    end
  end

  return ans
end



######## Generations of edge multiplicities

function multiedges_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, col_aut::Int64, max_genus::Int64, top_genus::Int64, Multi)
  
  return Iterators.flatmap(edge_mult_plus_edge_mult_aut -> all_multiedges_fixed_multi_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, edge_mult_plus_edge_mult_aut[1], edge_mult_plus_edge_mult_aut[2], max_genus, top_genus), unique_edge_multi_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, col_aut, Multi))
end

function all_multiedges_fixed_multi_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, edge_mult::Vector{Int64}, edge_mult_aut::Int64, max_genus::Int64, top_genus::Int64)
  
  genus_coming_from_multiedges = max_genus - top_genus - sum(gen_dist)
  # ### no need to implement multiedges:
  # (genus_coming_from_multiedges == 0) && return [([[edge_mult[i]] for i in eachindex(edge_mult)], edge_mult_aut)]
  
  return Iterators.flatmap(p -> multiedges_fixed_multi_and_part_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, edge_mult, edge_mult_aut, p), Combinatorics.partitions(genus_coming_from_multiedges + Graphs.ne(top_graph_Graphs), Graphs.ne(top_graph_Graphs)))
end

function multiedges_fixed_multi_and_part_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, edge_multi::Vector{Int64}, edge_multi_aut::Int64, p::Vector{Int64})

  ans = Set{Tuple{NTuple{Graphs.ne(top_graph_Graphs), Vector{Int64}}, Int64}}()

  edge_mult_dict = Dict{Edge, Int}(edges(top_graph) .=> edge_multi)
  
  for array_number_multiedges = Combinatorics.multiset_permutations(p, Graphs.ne(top_graph_Graphs)) # distribute the parts of p to the edges
    edge_mult_number_multiedges = Dict{Edge, Int}(edges(top_graph) .=> array_number_multiedges)
    
    any(e -> edge_mult_dict[e] < edge_mult_number_multiedges[e], edges(top_graph)) && continue # if any edge has more multiedges than its multiplicity, skip

    for multi_edges in Iterators.product([Combinatorics.partitions(edge_mult_dict[e], edge_mult_number_multiedges[e]) for e in edges(top_graph)]...)

      multi_edges_dict = Dict{Edge, Vector{Int64}}(edges(top_graph) .=> multi_edges)
      found = false
      
      for (multi_edges2, multiedges_aut) in ans

        if (edge_multi_aut > 1)

          multi_edges2_dict = Dict{Edge, Vector{Int64}}(edges(top_graph) .=> multi_edges2)
          
          vertex_rel(u, v) = (col[u] == col[v]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
          edge_rel(u, v) = multi_edges_dict[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == multi_edges2_dict[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
          
          if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel, edge_relation=edge_rel) # check isomorphism
            found = true
            # seen[color2[1]] -= 1 # use one copy of this coloring
            break
          end
        end
      end

      if !found # new edge_mult if not found
        
        aut = 1
        if edge_multi_aut > 1
          vertex_rel_2(u, v) = (col[u] == col[v]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
          edge_rel_2(u, v) = multi_edges_dict[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == multi_edges_dict[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
          aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel_2, edge_relation=edge_rel_2) # count automorphisms of the coloring
        end

        push!(ans, (multi_edges, aut))
      end

    end

  end

  return ans
end


function unique_edge_multi_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, gen_dist_aut, Multi)

  # Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(top_graph)], beta)

  gen_dist_aut == 1 && return Iterators.zip(Multi, Iterators.cycle([gen_dist_aut])) # if no coloring automorphisms, no need to check edge multiplicities

  ans = Set{Tuple{Vector{Int64}, Int64}}()
  
  for edgeMult_array in Multi # iterate edge multiplicities

    edgeMult = Dict{Edge, Int}(edges(top_graph) .=> edgeMult_array)
    found = false
    
    for (edgeMult_array2, edge_mult_aut) in ans

      edgeMult2 = Dict{Edge, Int}(edges(top_graph) .=> edgeMult_array2)
      
      vertex_rel(u, v) = (col[v] == col[u]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
      edge_rel(u, v) = edgeMult[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == edgeMult2[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
      
      if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel, edge_relation=edge_rel) # check isomorphism
        found = true
        # seen[color2] -= 1 # use one copy of this coloring
        break
      end
    end

    if !found # new edge_mult if not found
      
      aut = 1

      vertex_rel_2(u, v) = (col[v] == col[u]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
      edge_rel_2(u, v) = edgeMult[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == edgeMult[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
      aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel_2, edge_relation=edge_rel_2) # count automorphisms of the coloring


      push!(ans, (edgeMult_array, aut))
    end
  end

  return ans
  
end



######## Generations of marks

# function unique_marks_dist_fixed_WRC(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, m_inv::Vector{Int64}, multi_edges, multi_edges_aut::Int64, n_marks::Int64)

#   multi_edges_aut == 1 && return Iterators.zip(Combinatorics.multiset_permutations(m_inv, length(m_inv)), Iterators.cycle([multi_edges_aut])) # if no coloring automorphisms, no need to check edge multiplicities

#   ans = Set{Tuple{Vector{Int64}, Int64}}()

#   multi_edges_dict = Dict{Edge, Vector{Int64}}(edges(top_graph) .=> multi_edges)
#   edge_rel(u, v) = multi_edges_dict[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == multi_edges_dict[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check

#   for multiset_perm in Combinatorics.multiset_permutations(m_inv, n_marks)
    
#     found = false

#     for (multiset_perm2, multiset_perm_aut) in ans
      
#       # vertex_rel(u, v) = (col[u] == col[v]) && (gen_dist[u] == gen_dist[v]) && all(i -> (u == multiset_perm[i]) == (v == multiset_perm2[i]), 1:n_marks) # vertex relation for isomorphism check
#       vertex_rel(u, v) = (col[u] == col[v]) && (gen_dist[u] == gen_dist[v]) && all(i -> sort([i for i in 1:n_marks if u == multiset_perm[i]]) == sort([i for i in 1:n_marks if v == multiset_perm[i]]), 1:n_marks) # vertex relation for isomorphism check
      
#       if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel) # check isomorphism
#         found = true
#         # seen[color2] -= 1 # use one copy of this coloring
#         break
#       end
#     end

#     if !found # new coloring if not found
#       vertex_rel_2(u, v) = (col[u] == col[v]) && (gen_dist[u] == gen_dist[v]) && all(i -> (u == multiset_perm[i]) == (v == multiset_perm[i]), 1:n_marks) # vertex relation for isomorphism check
#       aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel_2, edge_relation=edge_rel) # count automorphisms of the coloring
#       push!(ans, (multiset_perm, aut))
#     end
#   end

#   return ans
# end

# ##### Generations of edge multiplicities with multiedges
# function unique_multi_edge_fixed_multi(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, col_aut::Int64, max_genus::Int64, top_genus::Int64, edge_mult_array::Vector{Int64})

#   ans = Set{Vector{Vector{Int64}}}()

#   edge_mult_dict = Dict{Edge, Int}(edges(top_graph) .=> edge_mult_array)
  
#   for array_number_multiedges = Combinatorics.with_replacement_combinations(p, Graphs.ne(top_graph)) # distribute the parts of p to the edges
    
#     edge_mult_number_multiedges = Dict{Edge, Int}(edges(top_graph) .=> array_number_multiedges)
    
#     any(e -> edge_mult_dict[e] < edge_mult_number_multiedges[e], edges(top_graph)) && continue # if any edge has more multiedges than its multiplicity, skip
    
#     multi_edges = [collect(Combinatorics.partitions(edge_mult_dict[e], edge_mult_number_multiedges[e])) for e in edges(top_graph)]

#     push!(ans, multi_edges)
#   end

#   return ans
  
# end

# function multi_edge_multi_mod_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, col_aut::Int64, max_genus::Int64, top_genus::Int64, Multi::Set{Vector{Int64}})

#   genus_coming_from_multiedges = max_genus - top_genus - sum(gen_dist)
#   ### no need to implement multiedges:
#   (genus_coming_from_multiedges == 0) && return map((edge_mult, edge_mult_aut)->([[edge_mult[i]] for i in eachindex(edge_mult)], edge_mult_aut), unique_edge_multi_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, col_aut, Multi))

#   all_multiplicities = collect(unique_edge_multi_mod_iso(top_graph_Graphs, top_graph, gen_dist, col, col_aut, Multi))
  
#   return Iterators.flatmap(p -> unique_multi_edge_fixed_part(all_multiplicities, p), Combinatorics.partitions(genus_coming_from_multiedges + Graphs.ne(top_graph_Graphs), Graphs.ne(top_graph_Graphs)))
# end

# function unique_multi_edge_fixed_part(top_graph_Graphs::Graphs.SimpleGraph{Int64}, top_graph::Graph{Undirected}, gen_dist::Vector{Int64}, col::Vector{Int64}, all_multiplicities, p)
  
#   ans = Set{Tuple{Vector{Vector{Int64}}, Int64}}()

#   for (array_mult, array_mult_aut) in all_multiplicities

#     edge_mult_dict = Dict{Edge, Int}(edges(top_graph) .=> array_mult)
    
#     for array_number_multiedges = Combinatorics.with_replacement_combinations(p, Graphs.ne(top_graph)) # distribute the parts of p to the edges
      
#       edge_mult_number_multiedges = Dict{Edge, Int}(edges(top_graph) .=> array_number_multiedges)
      
#       any(e -> edge_mult_dict[e] < edge_mult_number_multiedges[e], edges(top_graph)) && continue # if any edge has more multiedges than its multiplicity, skip
      
#       multi_edges = [collect(Combinatorics.partitions(edge_mult_dict[e], edge_mult_number_multiedges[e])) for e in edges(top_graph)]

#       found = false
      
#       for (array_mult2, aut) in ans

#         edge_mult_dict2 = Dict{Edge, Int}(edges(top_graph) .=> array_mult2)
#         vertex_rel_2(u, v) = (col[v] == col[u]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
#         edge_rel_2(u, v) = edge_mult_dict[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == edge_mult_dict2[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
        
#         if (array_mult_aut > 1) && Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel, edge_relation=edge_rel) # check isomorphism
#           found = true
#           # seen[color2] -= 1 # use one copy of this coloring
#           break
#         end
#       end

#       if !found # new edge_mult if not found
#         vertex_rel_2(u, v) = (col[v] == col[u]) && (gen_dist[v] == gen_dist[u]) # vertex relation for isomorphism check
#         edge_rel_2(u, v) = edge_mult_dict[Edge(max(Graphs.src(v), Graphs.dst(v)), min(Graphs.src(v), Graphs.dst(v)))] == edge_mult_dict[Edge(max(Graphs.src(u), Graphs.dst(u)), min(Graphs.src(u), Graphs.dst(u)))] # edge relation for isomorphism check
#         aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=vertex_rel_2, edge_relation=edge_rel_2) # count automorphisms of the coloring
#         push!(ans, (array_mult, aut))
#       end
#     end
#   end
# end