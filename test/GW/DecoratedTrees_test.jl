# example for P^3
P3 = GKMproj_space(3);

# build decorated tree, which is 1-2-3, mapping 1 -> 1, 2-> 4, 3->3, and marked points (1, 2) at 1, and (3) at 2.
tree = Graph{Undirected}(2);
add_edge!(tree, 1, 2);
#add_edge!(tree, 2, 3);
#add_edge!(tree, 1, 3);

#vDict = [1, 4, 3];
vDict = [1, 4]

# edgeMult = Dict(Edge(1, 2) => 1, Edge(2, 3) => 1);
#edgeMult = Dict(Edge(1, 2) => 1, Edge(1, 3) => 1);
edgeMult = Dict(Edge(1, 2) => 1)

marks = [1, 2]

dt = decoratedTree(P3, tree, vDict, edgeMult, marks);

P = ev(1, pointClass(1, P3)) * ev(2, pointClass(4, P3))

println("Contribution of Edge(1,4) with d=1 in P3:")
println(GWTreeContribution(dt, P))


println("Contribution of Edge(1,4) with d=2 in P3:")
edgeMult = Dict(Edge(1, 2) => 2)
dt = decoratedTree(P3, tree, vDict, edgeMult, marks);
println(GWTreeContribution(dt, P))


#Warning: the example below does not work because vertices(G) throws an error for graphs G without edges.
#Conclusion: We always need special treetment for beta = zero.
println("Contribution of single point as decorated tree:")
marks = [1, 1, 1]
tree = Graph{Undirected}(1)
vDict = [1]
P = class_one()
dt = decoratedTree(P3, tree, vDict, edgeMult, marks);
res = GWTreeContribution(dt, P)
println(res)