P3 = GKMproj_space(3);

c1 = pointClass(1, P3); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(1, P3); #gens(H)[1] # pointClass(1, R) # H([C(0), C(0), C(0), C(1)]) # gens(H)[2]
c3 = pointClass(4, P3); #gens(H)[4] # pointClass(4, R)

P = ev(1, c1)*ev(2, c3)*ev(3, c3);
println(integrateGKM(P3, edgeCurveClass(P3, Edge(1, 2)), 3, P))

P2 = GKMproj_space(2);

(t1, t2, t3) = gens(P2.equivariantCohomology.coeffRing);
c_1 = pointClass(1, P2);
c_2 = pointClass(2, P2);
c_3 = pointClass(3, P2);

P = ev(1, c_1)*ev(2, c_1)*ev(3, c_1)*ev(4, c_1)*ev(5, c_1)
println(integrateGKM(P2, 5*edgeCurveClass(P2, Edge(1, 2)), 5, P))

P = ev(1, c_1)*ev(2, c_1)
println(integrateGKM(P2, 2*edgeCurveClass(P2, Edge(1, 2)), 2, P))

#(-48*t1^6 + 144*t1^5*t2 + 144*t1^5*t3 - 80*t1^4*t2^2 - 560*t1^4*t2*t3 - 80*t1^4*t3^2 - 80*t1^3*t2^3 + 560*t1^3*t2^2*t3 + 560*t1^3*t2*t3^2 - 80*t1^3*t3^3 + 144*t1^2*t2^4 - 336*t1^2*t2^3*t3 - 336*t1^2*t2^2*t3^2 - 336*t1^2*t2*t3^3 + 144*t1^2*t3^4 - 80*t1*t2^5 + 112*t1*t2^4*t3 + 112*t1*t2^3*t3^2 + 112*t1*t2^2*t3^3 + 112*t1*t2*t3^4 - 80*t1*t3^5 + 16*t2^6 - 16*t2^5*t3 - 16*t2^4*t3^2 - 16*t2^3*t3^3 - 16*t2^2*t3^4 - 16*t2*t3^5 + 16*t3^6)//(t1^6 - 3*t1^5*t2 - 3*t1^5*t3 + t1^4*t2^2 + 13*t1^4*t2*t3 + t1^4*t3^2 + 3*t1^3*t2^3 - 13*t1^3*t2^2*t3 - 13*t1^3*t2*t3^2 + 3*t1^3*t3^3 - 2*t1^2*t2^4 - t1^2*t2^3*t3 + 21*t1^2*t2^2*t3^2 - t1^2*t2*t3^3 - 2*t1^2*t3^4 + 4*t1*t2^4*t3 - 7*t1*t2^3*t3^2 - 7*t1*t2^2*t3^3 + 4*t1*t2*t3^4 - 2*t2^4*t3^2 + 5*t2^3*t3^3 - 2*t2^2*t3^4)