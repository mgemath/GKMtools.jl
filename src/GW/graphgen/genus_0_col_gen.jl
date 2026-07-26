struct col_it
    ls::Vector{Int64}
    col_dict::Dict{Int64,Vector{Int64}}
    rev_dfs::Vector{Int64}
    rev_positions::Vector{Int64}
    parents::Vector{Int64}
    has_ci::Bool
    subgraph_ends_rev::Vector{Int64}
    subgraph_ends::Vector{Int64}
    left_siblings::Vector{Int64}
    next_colors::Dict{Tuple{Int64,Int64},Int64}
end

const Marks_type = Vector{Int64}
function col_it_init(ls::Vector{Int64}, col_dict::Dict{Int64,Vector{Int64}})

    _validate_level_sequence(ls)
    _validate_color_dictionary(col_dict)
    n::Int64 = length(ls)
    par = zeros(Int64, n)
    last_at_level = Dict{Int64,Int64}(ls[1] => 1)
    for v in 2:n
        par[v] = last_at_level[ls[v] - 1]
        last_at_level[ls[v]] = v
    end

    children = [Int64[] for _ in 1:n]
    for v in 2:n
        push!(children[par[v]], v)
    end

    stack = Int64[1]
    rev_dfs = Int64[]
    while !isempty(stack)
        v = pop!(stack)
        push!(rev_dfs, v)
        append!(stack, children[v])
    end
    rev_positions = invperm(rev_dfs)
    sub_end = _subtree_ends(ls)

    left_siblings = fill(Int64(-1), n)
    for siblings in children
        for i in 2:length(siblings)
            left_siblings[siblings[i]] = siblings[i - 1]
        end
    end

    next_colors = Dict{Tuple{Int64,Int64},Int64}()
    for (parent_color, allowed) in col_dict, i in eachindex(allowed)
        next_colors[(parent_color, allowed[i])] = i == length(allowed) ? 0 : allowed[i + 1]
    end

    return col_it(ls, col_dict, rev_dfs, rev_positions, par,
        my_has_central_involution(ls),
        _reverse_subtree_ends(ls, rev_dfs),
        sub_end, left_siblings, next_colors), par, sub_end
end

function _validate_color_dictionary(col_dict::Dict{Int64,Vector{Int64}})
    isempty(col_dict) && throw(ArgumentError("the color adjacency dictionary cannot be empty"))
    colors = Set(keys(col_dict))
    sort!(collect(colors)) == collect(1:length(colors)) ||
        throw(ArgumentError("colors must be numbered consecutively starting at 1"))
    for (color, neighbors) in col_dict
        isempty(neighbors) && throw(ArgumentError("color $color must have at least one allowed neighboring color"))
        issorted(neighbors) || throw(ArgumentError("allowed neighboring colors for color $color must be sorted"))
        length(unique(neighbors)) == length(neighbors) || throw(ArgumentError("allowed neighboring colors for color $color must be unique"))
        all(in(colors), neighbors) || throw(ArgumentError("color $color refers to a color missing from the dictionary"))
        color ∉ neighbors || throw(ArgumentError("loops in the color adjacency dictionary are not supported"))
        all(other -> color in col_dict[other], neighbors) || throw(ArgumentError("the color adjacency dictionary must be symmetric"))
    end
    return nothing
end

function Base.iterate(CI::col_it, col::Vector{Int64}=Int64[])

    if isempty(col)
        c::Vector{Int64} = first_coloring(CI.ls, CI.col_dict)
        return (c...,), c #otherwise (c...,), c
    end


    #=
        A vertex is colorable if it hasn't reached its maximum possible color.

        The vertex 1 is colorable if either 
        1.) The graph has no central involution and the color is 1 is not the max color; or 
        2.) The graph has central involutiuon and its color is not the second largest color. 
    =#

    # this variable shows how far we need to look in the tree to find the vertex whose 
    # color is modified
    look_end::Int64 = CI.subgraph_ends[1]

    # in colorable we keep the index of the vertex that is to be colored
    colorable::Int64 = 0 #colorable_1 ? 1 : 0

    # we go through the vertices in reverse DFS order
    k::Int64 = 2

    while k <= length(CI.ls)
        v::Int64 = CI.rev_dfs[k]

        # compute end point of the subgraph T_v stemming at vertex v
        v_end::Int64 = CI.subgraph_ends[v]
        v_end_rev::Int64 = CI.subgraph_ends_rev[v]
        # compute the root left sibling of the subgraph T_v
        # (the subgraph stemming from a vertex with the same parent as v, left of v)
        left_s::Int64 = CI.left_siblings[v]

        left_sibling::UnitRange{Int64} = 1:0

        if left_s != -1
            # the left sibling exists
            # we compute the range of indices that correspond to the left sibling
            left_sibling = left_s:(v-1)
        end

        # check if T_v is isomorphic to the left sibling and if it has the same coloring
        # If yes, then the coloring of T_v cannot be increased and hence this part of 
        # the tree can be ignored

        if view(CI.ls, v:v_end) == view(CI.ls, left_sibling) && view(col, v:v_end) == view(col, left_sibling)
            k = CI.rev_positions[left_s]
            continue
        end

        # we decide if the color of v can be increased

        # count_if_increaseble = count(i -> i > col[v], CI.col_dict[col[CI.parents[v]]])
        if col[v] < CI.col_dict[col[CI.parents[v]]][end]
            # if (col[v] < CI.max_col - 1) || (col[v] == CI.max_col - 1 && col[CI.parents[v]] != CI.max_col)
            # we update colorable if necessary
            if colorable < v
                colorable = v
            end
            #we update look_end
            look_end = v_end_rev
        end

        if v == look_end && colorable != 0
            # we got the end of a branch that has colorable vertex
            # we need to look no further
            break
        end
        k += 1
    end

    if colorable == 0 # if all vertices >1 are not colorable
        new_root_color = 0  # if necessary, we need to increase the color of ls[1]
        if CI.has_ci # if g has the bad automorphism
            for c in (col[1]+1):maximum(keys(CI.col_dict))#CI.max_col # find the first color, greater than col[1], which is not maximal among its companion colors
                if c < CI.col_dict[c][end]
                    new_root_color = c
                    break
                end
            end
        else
            if col[1] < maximum(keys(CI.col_dict))#CI.max_col  # if g has not the bad automorphism, we can increase col[1] safely by 1
                new_root_color = col[1] + 1
            end
        end

        if new_root_color == 0 # none of the above conditions is satisfied, then we cannot incrase the coloration
            return nothing
        else
            col[1] = new_root_color  # we increase the color of ls[1]
            colorable = 1
        end
    else
        col[colorable] = next_color(CI, col, colorable)
    end


    # we reset the coloring of each subtree on the right side of v
    # starting from v+1

    j::Int64 = colorable + 1
    while j <= length(CI.ls)
        color_j = (CI.has_ci && j == 2) ? CI.col_dict[col[1]][findfirst(i -> i > col[1], CI.col_dict[col[1]])] : CI.col_dict[col[CI.parents[j]]][1]

        # find the end of the subtreee T_j stemming from vertex j
        es::Int64 = CI.subgraph_ends[j]

        # compute the minimal coloring for the subtree T_j
        # with parent color being the color of parent[j]
        col[j:es] = first_coloring(CI.ls[j:es], CI.col_dict, color_j; validate=false)
        # col[j:es] = my_minimal_coloring2(CI.ls[j:es], CI.col_dict, parent_color=col[CI.parents[j]], root_color=root_color)
        # the following subtree that needs to be dealt with stems from the end of T_j plus 1 
        j = es + 1
    end

    # return the coloring computed
    return (col...,), col

end

function first_coloring(ls::Vector{Int64}, col_dict::Dict{Int64,Vector{Int64}}, color_root::Int64=1; validate::Bool=true)::Vector{Int64}

    if validate
        _validate_level_sequence(ls; require_root_level=false)
        _validate_color_dictionary(col_dict)
    end
    haskey(col_dict, color_root) || throw(ArgumentError("root color $color_root is missing from the color dictionary"))
    ans::Vector{Int64} = fill(color_root, length(ls))
    last_at_level = Dict{Int64,Int64}(ls[1] => 1)
    for i in 2:length(ls)
        parent = last_at_level[ls[i] - 1]
        ans[i] = col_dict[ans[parent]][1]
        last_at_level[ls[i]] = i
    end

    return ans
end

function next_color(CI::col_it, col::Vector{Int64}, v::Int64)
    return CI.next_colors[(col[CI.parents[v]], col[v])]
end


#########AUX FUNCTIONS#######

function my_root_of_left_sibling_subtree(par::Vector{Int64}, v::Int64)::Int64

    p::Int64 = par[v]
    if p == 0
        return -1
    end
    desc::Vector{Int64} = [x for x in my_children_vertices(par, p) if x < v]

    if isempty(desc)
        return -1
    end

    return maximum(desc)
end

function my_children_vertices(par::Vector{Int64}, v::Int64)::Vector{Int64}

    return findall(i -> par[i] == v, eachindex(par))
end

function _reverse_subtree_ends(ls::Vector{Int64}, order::Vector{Int64})
    n = length(order)
    result = Vector{Int64}(undef, n)
    open = Int64[]
    for pos in n:-1:1
        level = ls[order[pos]]
        while !isempty(open) && ls[order[open[end]]] > level
            pop!(open)
        end
        stop = isempty(open) ? n : open[end] - 1
        result[order[pos]] = order[stop]
        push!(open, pos)
    end
    return result
end

function _subtree_ends(ls::Vector{Int64})
    ends = fill(length(ls), length(ls))
    open = Int64[]
    for v in eachindex(ls)
        while !isempty(open) && ls[v] <= ls[open[end]]
            ends[pop!(open)] = v - 1
        end
        push!(open, v)
    end
    return ends
end

function my_end_of_subgraph(ls::Vector{Int64}, v::Int64)

    # find the first vertex whose level is not grater then the level of v
    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return findfirst(i -> i == length(ls) || (i >= v && ls[i+1] <= ls[v]), eachindex(ls))
end

function my_end_of_subgraph_rev(ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64)::Int64

    # position of v in r_dfs 
    ps = findfirst(x -> r_dfs[x] == v, 1:length(ls))
    end_ver::Union{Nothing,Int64} = findfirst(k -> ls[k] <= ls[v], view(r_dfs, (ps+1):length(r_dfs)))

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end

function end_of_subgraph_rev(ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64)::Int64 #to review

    # position of v in r_dfs 
    ps = findfirst(x -> r_dfs[x] == v, 1:length(ls))
    end_ver::Union{Nothing,Int64} = findfirst(k -> ls[k] <= ls[v], r_dfs[(ps+1):end])

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end

function my_has_central_involution(ls::Vector{Int64})::Bool

    # if the graph is o--o then the answer is yes
    if ls == [1, 2]
        return true
    end

    if iseven(length(ls))
        two = findfirst(i -> i > 2 && ls[i] == 2, eachindex(ls))
        if two !== nothing && view(ls, 3:(two-1)) == 1 .+ view(ls, two:length(ls))
            return true
        end
    end

    return false
end

######Of Iterators#########

Base.eltype(::Type{col_it}) = Tuple{Vararg{Int64}}
Base.IteratorSize(::Type{col_it}) = Base.SizeUnknown()

### this function counts the isomorphisms of a tree with level sequence ls. Optionally, the tree can be colored with coloration col

function count_iso(ls::Vector{Int64}, col::Tuple{Vararg{Int64}}, marks::Marks_type)::Int64

    isempty(marks) && return count_iso(ls, col)

    temp_col = Vector{Int64}(undef, length(ls))
    for i in eachindex(temp_col)
        if i in marks
            temp_col[i] = typemax(Int64) - i
        else
            temp_col[i] = col[i]
        end
    end

    return count_iso(ls, (temp_col...,))
end

function count_iso(ls::Vector{Int64}, col::Tuple{Vararg{Int64}}=())::Int64

    is_empty::Bool = isempty(col)

    if length(ls) < 3
        if is_empty
            return ls[1] == 1 ? length(ls) : 1 # ls[1] == 1 is equivalent to be a starting-point graph
        else
            return 1
        end
    end

    my_child = findall(i -> i > length(ls) || ls[i] == ls[1] + 1, 1:(length(ls)+1))  # find all the children of the root, plus length(ls)+1


    if is_empty
        if ls[1] == 1 # ls[1] == 1 is equivalent to be a starting-point graph
            if iseven(length(ls))  # necessary condition to have the bad involution
                if view(ls, (my_child[1]+1):(my_child[2]-1)) == 1 .+ view(ls, my_child[2]:length(ls))  # check if it has the bad involution
                    a = count_iso(ls[my_child[1]:(my_child[2]-1)])   # compute the number of colorations of the main subgraph
                    return 2 * (a^2)
                end
            end
        end
    end

    # here we compute the number of colorations of each subgraph
    last_sub::Vector{Int64} = ls[my_child[1]:(my_child[2]-1)]  # this is the subgraph
    last_col::Tuple{Vararg{Int64}} = col

    if !is_empty
        last_col = col[my_child[1]:(my_child[2]-1)]
    end

    m::Int64 = 1
    ans::Int64 = 1


    for x in 3:length(my_child)
        if last_sub == view(ls, my_child[x-1]:(my_child[x]-1)) && (is_empty || last_col == col[my_child[x-1]:(my_child[x]-1)]) #view(col, my_child[x-1]:my_child[x]-1))
            m += 1
        else
            a = count_iso(last_sub, last_col)   # number of colorations of this subgraph...
            ans *= (a^m) * factorial(m)
            last_sub = ls[my_child[x-1]:(my_child[x]-1)]  # pass to the next subgraph
            # last_col = col[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
            if !is_empty
                last_col = col[my_child[x-1]:(my_child[x]-1)]
            end
            m = 1
        end
    end

    a = count_iso(last_sub, last_col)  # compute the number of colorations of the last subgraph
    ans *= (a^m) * factorial(m)       # with multiplicity

    return ans
end

#This functions computes the number of colorations. It is not called anywhere in the package, but I decided to keep it because it could be useful in the future

# function n_colorations(ls::Vector{Int64}, n::Int64, magn::Int64)::Int64  # magn is the number of neighbors cones of a cone. It is nc[1]

#     if length(ls) < 3
#         if ls[1] == 1
#             if length(ls) == 1
#                 return n
#             else
#                 return div(n*magn, 2)
#             end
#         else
#             return magn^length(ls)
#         end
#     end

#     my_child = findall(i -> i>length(ls) || ls[i]==ls[1]+ 1, 1:length(ls)+1)  # find all the children of the root, plus length(ls)+1


#     if ls[1] == 1  # ls[1] == 1 is equivalent to be a starting-point graph
#         if iseven(length(ls))  # necessary condition to have the bad involution
#             if view(ls,my_child[1]+1:my_child[2]-1) == 1 .+ view(ls,my_child[2]:length(ls))  # check if it has the bad involution
#                 c = n_colorations(ls[my_child[1]:my_child[2]-1], n, magn)   # compute the number of colorations of the main subgraph
#                 return div(n*(c^2),2*magn)  # return the number of coloration computing only the number of colorations of the main subtree
#             end
#         end
#         ans = n  # in the starting-point graph, the root can assume n values since it has no parent
#     else
#         ans = magn # we are not in the starting-point graph, so the root has a parent to deal with
#     end

#     # here we compute the number of colorations of each subgraph
#     last_sub = ls[my_child[1]:my_child[2]-1]  # this is the subgraph
#     m = 1                                     # this is its multiplicity

#     for x in 3:length(my_child)
#         if last_sub == view(ls, my_child[x-1]:my_child[x]-1)  # we run up to find a different subgraph
#             m += 1
#         else
#             c = n_colorations(last_sub, n, magn)   # number of colorations of this subgraph...
#             ans *= binomial(c+m-1,m)         # ...counted with multiplicity
#             last_sub = ls[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
#             m = 1
#         end
#     end

#     c = n_colorations(last_sub, n, magn)  # compute the number of colorations of the last subgraph
#     ans *= binomial(c+m-1,m)        # with multiplicity

#     return ans
# end