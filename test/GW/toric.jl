P3 = projective_space(NormalToricVariety,3);#thespaceP^3.
X = gkm_graph(domain(blow_up(P3,ideal(gens(cox_ring(P3))[1:3]))));

for e in edges(X.g)
    println("$(X.labels[src(e)]) -> $(X.labels[dst(e)]) => ", edgeCurveClass(X, e))
end

# (-1, 1)
edgeCurveClass(X, edgeFromLabels(X, "1", "2"))
edgeCurveClass(X, edgeFromLabels(X, "1", "3"))
edgeCurveClass(X, edgeFromLabels(X, "2", "3"))

#(1, 0)
edgeCurveClass(X, edgeFromLabels(X, "1", "4"))


H=edgeCurveClass(X, edgeFromLabels(X, "5", "4")); #(0,1)
# E=edgeCurveClass(H2, edgeFromLabels(X, "1", "4")) #(1,0)
E=edgeCurveClass(X, edgeFromLabels(X, "1", "2")) #(-1,1)
d=2;#thiscanbeanynonnegativeinteger
e=-1;#thiscanbeanynonpositiveinteger
beta=d*H+e*E;
# beta = 3*E


P= prod(i -> ev(i, pointClass(1, X)), 1:(2*d+e));
println(integrateGKM(X, beta, 2*d+e, P));
