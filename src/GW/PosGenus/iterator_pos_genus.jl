function geng(g::Int64, n::Int64)::Base.EachLine{IOBuffer}

    n < g + 2 && return eachline(IOBuffer("")) # skip impossible cases
    n_edge = g + n - 1 # number of edges for connected graph with genus g and n vertices
    cmd = nauty_jll.geng_path * " -c -q $n $n_edge:$n_edge"
    iter_of_g6 = eachline(IOBuffer(read(`sh -c $cmd`, String)))

    return iter_of_g6
end

function compute_aut(M::Matrix{Int64})::Int64

    # graph_of_Graphs = Graphs.SimpleGraph{Int64}(M)

    # return Graphs.Experimental.count_isomorph(graph_of_Graphs, graph_of_Graphs)
    return 1 # TO BE FIXED LATER
end

function graph6_to_adjacency_matrix(s::String)::Matrix{Int64}
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
    adj = zeros(Int, n, n)

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
                    adj[i_val+1, j_val+1] = 1
                    adj[j_val+1, i_val+1] = 1
                end
                bit_index += 1
            end
        end
    end

    return adj
end
