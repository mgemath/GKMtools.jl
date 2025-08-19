using Test, Combinatorics, Oscar, GKMtools

#some working examples
M = free_module(ZZ, 2);
g = Graph{Undirected}(3);
add_edge!(g, 1, 2);
add_edge!(g, 2, 3);
add_edge!(g, 1, 3);
de = Dict(edges(g) .=> [zero(M) for _ in n_edges(g)]);
labs = ["a", "b", "c"];
G = gkm_graph(g, labs, M, de)
[G, G]

de2 = Dict(edges(g) .=> [gens(M)[1] for _ in n_edges(g)]);
G2 = gkm_graph(g, labs, M, de2)

GKMproj_space(4)
GKMproj_space(4, label = "a")

@req is3_indep(GKMproj_space(4)) "Fail1"

#Figure 1.(b)
labs2 = String[]

for (i,j,k) in Iterators.product(1:3, 1:3, 1:3)
    (i == j || j == k || i == k) && continue
    "$i$j$k" in labs2 && continue
    push!(labs2, "$i$j$k")
end
G2 = empty_gkm_graph(6, 3, labs2)
N = G2.M;
x1 = gens(N)[1];
x2 = gens(N)[2];
x3 = gens(N)[3];

GKMadd_edge!(G2, "123", "213", x1-x2);
GKMadd_edge!(G2, "123", "321", x1-x3);
GKMadd_edge!(G2, "123", "132", x2-x3);

GKMadd_edge!(G2, "132", "231", x1-x2);
GKMadd_edge!(G2, "132", "312", x1-x3);

GKMadd_edge!(G2, "312", "321", x1-x2);
GKMadd_edge!(G2, "312", "213", x2-x3);

GKMadd_edge!(G2, "321", "213", x2-x3);

GKMadd_edge!(G2, "213", "213", x1-x3);

is3_indep(G2)

F=flag_gkm_graph([1,1,1]) # this is the same as the variety G2, or Figure 1.(b). It is the variety of complete flags in C^3

flag_gkm_graph([1, 3]) # this is the projective space of dimension 3-1

FlagVar34 = flag_gkm_graph([3, 4]) # this is the Grassmannian G(k,n) where k=3 and n = 3+4
println("(Combinatorial) Betti numbers of flag variety with quotient dimensions (3, 4): $(bettiNumbers(FlagVar34))")
@req GKM_isValid(FlagVar34) "Fail2"

h=hirzebruch_surface(NormalToricVariety, 6);

Gh = gkm_graph(h) # from a toric variety to a GKM graph
@req bettiNumbers(Gh) == [1, 2, 1] "Fail3"
@req GKM_isValid(Gh) "Fail3"

include("GKMsubgraph_test.jl")