# Case study of the GKM graph of the full flag variety F3 for C^3 but with different connection

R = root_system(:A, 2);
F3_alt_con = generalized_gkm_flag(R);
#S = QH_structure_constants(F3_alt_con);
con = get_connection(F3_alt_con);

println(con)
#hDict1 = [GKMtools._h(Edge(1, 2), d, con, F3_alt_con.equivariantCohomology) for d in 1:20];

# connection for (1,2)
println(con.a[Edge(1, 2), Edge(1, 4)])
con.con[Edge(1, 2), Edge(1, 4)] = Edge(2, 6)
con.con[Edge(2, 1), Edge(2, 6)] = Edge(1, 4)
con.con[Edge(1, 2), Edge(1, 5)] = Edge(2, 3)
con.con[Edge(2, 1), Edge(2, 3)] = Edge(1, 5)
# connection for (1 5)
con.con[Edge(1, 5), Edge(1, 4)] = Edge(5, 3)
con.con[Edge(5, 1), Edge(5, 3)] = Edge(1, 4)
con.con[Edge(1, 5), Edge(1, 2)] = Edge(5, 6)
con.con[Edge(5, 1), Edge(5, 6)] = Edge(1, 2)
conNew = build_GKM_connection(F3_alt_con, con.con)
set_connection!(F3_alt_con, conNew)

con = get_connection(F3_alt_con)
println(con.a[Edge(1, 2), Edge(1, 4)])

#hDict2 = [GKMtools._h(Edge(1, 2), d, con, F3_alt_con.equivariantCohomology) for d in 1:20];


#BO = get_bruhat_order_of_generalized_flag(R);
#S = schubert_classes(F3, BO);
#println("QH structure constants of F3 with alternate connection:")
#S2 = QH_structure_constants(F3_alt_con; refresh=true)
#NOTE 1: The resulting QH is polynomial, homogeneous, associative, and commutative.
#NOTE 2: We get S == S2 is true! 

#SB = QH_structure_constants_in_basis(F3, S)
#schubertClassLabels = ["\\sigma_{" * BO.labels[i] * "}" for i in 1:length(BO.labels)];
#schubertClassLabels = ["\\sigma_{id}", "\\sigma_{s_1}", "\\sigma_{s_2s_1}", "", "\\sigma_{s_2}", "\\sigma_{s_1s_2}"]
#QH_print_structure_constants_in_basis(F3, S, schubertClassLabels)

# after the output, remove lines which are multiplication by one.