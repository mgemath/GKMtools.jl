
function convex_not_signed_example()
    G = empty_gkm_graph(8, 2, ["$i" for i in 1:8])
    x, y = gens(G.M)

    # add edges with weight x
    for (i, j) in [(1, 2), (4, 3), (5, 6), (8, 7)]
        add_edge!(G, i, j, x)
    end

    # add edges with weight y
    for (i, j) in [(1, 8), (2, 3), (7, 6), (4, 5)]
        add_edge!(G, i, j, y)
    end

    # add remaining edges
    add_edge!(G, 1, 4, x+2*y)
    add_edge!(G, 8, 5, x+2*y)
    add_edge!(G, 2, 7, -x+2*y)
    add_edge!(G, 3, 6, -x+2*y)

    return G
end