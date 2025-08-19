G = flag_gkm_graph([1, 3]);

S = GKMsubgraph_from_vertices(G, [1]);
blowupSub = blowupGKM(S);

X = blowupSub.super;
H = edgeCurveClass(X, edgeFromLabels(X, "2", "3"));
E = edgeCurveClass(X, edgeFromLabels(X, "[1>3]", "[1>2]"));

d=2;#thiscanbeanynonnegativeinteger
e=-1;#thiscanbeanynonpositiveinteger
beta=d*H+e*E;
# beta = 3*E

P= prod(i -> ev(i, pointClass(1, X)), 1:(2*d+e));
println(integrateGKM(X, beta, 2*d+e, P, show_bar = false));
