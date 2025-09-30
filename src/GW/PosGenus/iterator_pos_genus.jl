# const libsizeaut = joinpath(@__DIR__, "../../..", "deps", "libsizeaut." * Base.Libc.Libdl.dlext)
# function compute_aut(M::Matrix{Cint}, n::Int64)::Culong

#   return ccall((:sizeaut, libsizeaut), Culong, (Ptr{Cint}, Cint), M, n)
# end

function compute_aut(top_graph_Graphs::Graphs.SimpleGraph{Int64})::Int64
  return Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs)
end

function my_geng(g::Int64, n::Int64)::Base.EachLine{IOBuffer}

  # (g > 0) && n < ceil(Int, (3 + sqrt(1 + 8*g)) / 2) && return eachline(IOBuffer("")) # skip impossible cases, already checked in main loop
  n_edge = g + n - 1 # number of edges for connected graph with genus g and n vertices
  cmd = nauty_jll.geng_path * " -c -q $n $n_edge:$n_edge"
#   iter_of_g6 = eachline(IOBuffer(read(`sh -c $cmd`, String)))
  iter_of_g6 = eachline(IOBuffer(read(`sh -c $cmd`)))

  return iter_of_g6
end

function graph6_to_adjacency_matrix(s::String)::Matrix{Cint}
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
  adj = zeros(Cint, n, n)

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
          adj[i_val+1, j_val+1] = Cint(1)
          adj[j_val+1, i_val+1] = Cint(1)
        end
        bit_index += 1
      end
    end
  end

  return adj
end

function colorings_modulo_iso(top_graph_Graphs::Graphs.SimpleGraph{Int64}, nc::Dict{Int64,Vector{Int64}}, top_aut::Int64)::Base.Iterators.Flatten{Vector{Set{Tuple{Vector{Int64}, Int64}}}}
  return Iterators.flatten([unique_col_fixed_combination_w_counting(top_graph_Graphs, nc, top_aut, comb) for comb in Combinatorics.with_replacement_combinations(1:length(nc), Graphs.nv(top_graph_Graphs))])   
end

function unique_col_fixed_combination_w_counting(top_graph_Graphs::Graphs.SimpleGraph{Int64}, nc::Dict{Int64,Vector{Int64}}, top_aut::Int64, comb::Vector{Int64})::Set{Tuple{Vector{Int64}, Int64}}
  ans = Set{Tuple{Vector{Int64}, Int64}}()
  seen = Dict{Vector{Int64}, Int64}()

  for c in Combinatorics.multiset_permutations(comb, length(comb))

    all(e -> c[Graphs.src(e)] in nc[c[Graphs.dst(e)]], Graphs.edges(top_graph_Graphs)) || continue # check if coloring is valid
    
    found = false
    
    for color2 in ans

      seen[color2[1]] == 0 && continue # already used all the copies of this coloring
      color_rel(u, v) = (c[u] == color2[1][v]) # vertex relation for isomorphism check
      
      if Graphs.Experimental.has_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_rel) # check isomorphism
        found = true
        seen[color2[1]] -= 1 # use one copy of this coloring
        break
      end
    end

    if !found # new coloring if not found
      color_rel_2(u, v) = (c[u] == c[v])
      aut = Graphs.Experimental.count_isomorph(top_graph_Graphs, top_graph_Graphs, vertex_relation=color_rel_2) # count automorphisms of the coloring
      push!(ans, (c, aut))
      seen[c] = div(top_aut, aut) - 1 # how many copies of this coloring are there
    end

  end
  return ans
end

function genus_distribution(max_genus::Int64, genus::Int64, n_vert::Int64)::Vector{Vector{Int64}}

  max_genus == genus && return [zeros(Int64, n_vert)] # no distributions possible
  
  ans = Vector{Vector{Int64}}()

  for p in Combinatorics.partitions(max_genus - genus) 
    extended_p = vcat(p, zeros(Int64, n_vert - length(p)))
    for perm_extended_p in Combinatorics.multiset_permutations(extended_p, n_vert)
      push!(ans, perm_extended_p)
    end
  end

  return ans
end