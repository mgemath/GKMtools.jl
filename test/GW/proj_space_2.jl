# This is dedicated to computations using P2

P2 = GKMproj_space(2);
d = 3; # keep this <4, otherwise running time is huge
beta = d*edgeCurveClass(P2, Edge(1,2))

P = prod(i -> ev(i, pointClass((i % 3) + 1, P2)), 1:(3*d-1));
println(integrateGKM(P2, beta, 3*d-1, P))


# P = ev(1, c1)*ev(2, c1);
# @test integrateGKM(P2, 2, P, 1) == 1

# P = ev(1, c1)*ev(2, c1)*ev(3, c1)*ev(4, c1)*ev(5, c1);
# @test integrateGKM(P2, 5, P, 2) == 1

# P = class_one();
# @test integrateGKM(P2, 2, P, 1) == 0
# @test integrateGKM(P2, 0, P, 1) == 0

# P = ev(1, c1)*ev(2, c1)*ev(3, c1)*ev(4, c1)*ev(5, c1)*ev(6, c1)*ev(7, c1)*ev(8, c1);
# @test integrateGKM(P2, 8, P, 3) == 12