F = flag_gkm_graph([2,2]);
c1 = pointClass(1, F);

# P = class_one();
# @test integrateGKM(F, 0, P, 1) == 0 #expected 0

# P = class_one();
# @test integrateGKM(F, 1, P, 1) == 0 #expected 0

P = ev(1, c1)*Psi(2);
H2 = GKM_second_homology(F)
x=integrateGKM(F, 1*gens(H2.H2)[1], 1, P) #expected ?

println(x)

P = ev(1, c1)*Psi(2, 1);
x=integrateGKM(F, 1*gens(H2.H2)[1], 1, P) #expected ?
