# Hirzebruch surface with parameter a in ZZ.
a = 1 # a = 1 and 0 are the only positive cases.

H = gkm_graph_of_toric(hirzebruch_surface(NormalToricVariety, a))
print_curve_classes(H)
QH_structure_constants(H)

B = [poincare_dual(gkm_subgraph_from_edges(H, [e])) for e in [Edge(1, 2), Edge(2, 3), Edge(3, 4), Edge(4, 1)]]
M = Array{Any}(undef, (4, 4))
for i in 1:4, j in 1:4
    M[i, j] = B[i][j]
end
e = [QH_class(H, M[i,:]) for i in 1:4]
lab = ["x12, x23, x34, x41"]
