############################
# Betti numbers
############################

# @doc raw"""
#     betti_number(G::AbstractGKM_graph, i::Int) -> ZZRingElem

# Compute the `i`-th combinatorial Betti number of the GKM graph `G`, 
# as defined in [Guillemin--Zara, section 1.3](@cite). The combinatorial 
# Betti numbers equal the Betti numbers of the underlying GKM space if the 
# torus action is Hamiltonian.
# # Examples
# ```jldoctest
# julia> G = grassmannian(GKM_graph, 1, 3);

# julia> betti_number(P3,0)
# 1

# julia> betti_number(P3, 1)
# 0
# ```
# """
# function betti_number(G::AbstractGKM_graph, i::Int)
  
#   # check input
#   d = valency(G)::Int

#   if !(0 <= i <= 2*d) || isodd(i)
#       return ZZRingElem(0)
#   end

#   return betti_numbers(G)[(i >> 1) + 1]

# end
##############################################################################
# @doc raw"""
#     betti_numbers(G::AbstractGKM_graph) -> Vector{ZZRingElem}

# Compute all the even combinatorial Betti numbers of the GKM graph `G` as an array.
# # Examples
# ```jldoctest
# julia> P3 = projective_space(NormalToricVariety, 3)
# Normal toric variety

# julia> betti_number(P3,0)
# 1

# julia> betti_number(P3, 1)
# 0
# ```
# """
# function betti_numbers(G::AbstractGKM_graph)

#   for counter in 1:10^8 # arbitrary maximum number of attempts to avoid "while true"
  
#   xi = ZZ.(rand(Int, rank_torus(G)))  #TODO: find something without using random numbers
#     wxi = Dict{Edge, ZZRingElem}() # wxi stands for weight[e](xi)
#     isPolarizing = true
  
#     # calculate weight[e](xi) for all edges e
#     for e in edges(G.g)
  
#       wxi[e] = wxi[reverse(e)] = 0
  
#       for j in 1:rank_torus(G)
#         wxi[e] += xi[j] * G.w[e][j]
#         wxi[reverse(e)] += xi[j] * G.w[reverse(e)][j]
#       end
  
#       if wxi[e] == 0 || wxi[reverse(e)] == 0
#         isPolarizing = false
#         break
#       end
#     end
#     (!isPolarizing) && continue
  
#     # from here on, xi is known to be polarizing.
#     #betti = Dict{Int, Int}([i => 0 for i in 0:valency(G)])
#     betti = zeros(ZZRingElem, valency(G)+1)
  
#     for v in 1:n_vertices(G.g)
  
#       i = count(w -> wxi[Edge(v,w)] < 0, all_neighbors(G.g, v) )
#       betti[i+1] += ZZ(1)
#     end
#     return betti
#   end
# end