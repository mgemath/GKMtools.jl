
function get_vertex_inj_example()
    G = empty_gkm_graph(8, 3, ["v", "a1", "a2", "a3", "b1", "b2", "b3", "w"])
    e = gens(G.M)
    
    for i in 1:3
        add_edge!(G, "v", "a$i", e[i])
        add_edge!(G, "w", "b$i", e[i])
        for j in 1:3
            i == j && continue
            add_edge!(G, "a$i", "b$j", e[i] - e[j])
        end
    end

    return G
end