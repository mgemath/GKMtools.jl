# A2 = root_system(:A, 2);
# # generalized_flag(A2)
# SA2 = [simple_roots(A2)[1]];
# d = generalized_flag(A2, SA2)

# generalized_flag(root_system(:D, 4))
# generalized_flag(root_system(:G, 2))
# generalized_flag(root_system(:F, 4))
# C3 = root_system(:C, 3);
# SC31 = [simple_roots(C3)[3]];
# generalized_flag(C3, [simple_roots(C3)[3]])

# check if I am using the same as Oscar notation

# _generator_matrix(C3);
# _generator_matrix(A2);
# X=_generator_matrix(root_system(:A, 4))
# _generator_matrix(root_system(:B, 3))
# _generator_matrix(root_system(:C, 3))
# _generator_matrix(root_system(:D, 4))
# _generator_matrix(root_system(:F, 4))
# _generator_matrix(root_system(:G, 2))
# _generator_matrix(root_system(:E, 8))
# _generator_matrix(root_system(:E, 7))
# _generator_matrix(root_system(:E, 6))
# generalized_flag(root_system(:G, 2))
# rf = root_system(:F, 4)
# generalized_flag(rf, simple_roots(rf)[1:2])

# r6 = root_system(:E, 6)
# generalized_flag(r6, simple_roots(rf)[1:4])

# R = root_system([(:A, 2), (:A, 2)])
# generalized_flag(R)

# mapreduce(length, +, [[1,2],[1,2]]) # == 1 + 4 + 9

# ncC3 = root_system(matrix(ZZ, [2 -1 -2; -1 2 0; -1 0 2]))
# root_system(cartan_type_with_ordering(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2])) #([(:B, 2), (:C, 2)], [1, 3, 2, 4])
# nc = root_system(cartan_type(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2]))
# generalized_flag(nc)
# ([(:B, 2), (:C, 2)], [1, 3, 2, 4]))
# NC=sub(cartan_matrix(root_system([(:A, 2), (:A, 2)])), [1,3,2,4] , [1,3,2,4]);

# generalized_flag(root_system(NC));

# NC=sub(cartan_matrix(root_system([(:A, 1), (:A, 1)])), [2,1] , [2,1]);
# generalized_flag(root_system(NC))

# generalized_flag(root_system([(:A, 1), (:A, 2)]))
NC=sub(cartan_matrix(root_system([(:A, 1), (:A, 2)])), [3,1,2] , [3,1,2]);
generalized_flag(root_system(NC))

# generalized_flag(root_system([(:A, 1), (:A, 4)])); #0.510605 seconds (4.39 M allocations: 196.455 MiB, 13.31% gc time)
# generalized_flag2(root_system([(:A, 1), (:A, 4)])); #0.773291 seconds (6.86 M allocations: 333.550 MiB, 17.70% gc time)

# generalized_flag(root_system([(:A, 1)]))
# generalized_flag2(root_system([(:A, 1)]))
