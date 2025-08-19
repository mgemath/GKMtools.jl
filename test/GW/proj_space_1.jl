# This is dedicated to computations using P1

P1 = GKMproj_space(1);
(t1, t2) = gens(P1.equivariantCohomology.coeffRing);

c1 = pointClass(1, P1); #gens(H)[1] # pointClass(1, R) # H([t1, t2, t3, t4])
c2 = pointClass(2, P2);
beta = edgeCurveClass(P1, Edge(1, 2))

P = class_one();
@test integrateGKM(P1, beta, 0, P) == 1 #expected 1

P = ev(1, c1);
@test integrateGKM(P1, beta, 1, P) == 1

P = ev(1, c1)*ev(2, c1)*ev(3, c2)
@test integrateGKM(P1, beta, 3, P) == 1

@test integrateGKM(P1, 3*beta, 0, class_one()) == 0
@test integrateGKM(P1, 3*beta, 2, class_one()) == 0

P = Psi(1)
@test integrateGKM(P1, beta, 1, P) == -2

P = Psi(1, 1)
@test integrateGKM(P1, 1*beta, 1, P) == 2
P = Psi(2, 2)
@test integrateGKM(P1, 2*beta, 1, P) == 5//4