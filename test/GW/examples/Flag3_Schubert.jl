# Case study of the full flag variety F3 for C^3

R = root_system(:A, 2);
F3_schubert = generalized_gkm_schubert(R, "s1*s2");
BO = get_bruhat_order_of_generalized_flag(R);
S = schubert_classes(F3_schubert, BO);
QH_structure_constants(F3_schubert.self);
SB = QH_structure_constants_in_basis(F3_schubert.self, S)
schubertClassLabels = ["\\sigma_{id}", "\\sigma_{s_1}", "\\sigma_{s_2}", ""]
println("QH structure constants for Z_{s_1*s_2}:")
QH_print_structure_constants_in_basis(F3_schubert.self, S, schubertClassLabels)

F3_schubert_small = generalized_gkm_schubert(R, "s1");
S_small = schubert_classes(F3_schubert_small, BO);
QH_structure_constants(F3_schubert_small.self);
SB = QH_structure_constants_in_basis(F3_schubert_small.self, S_small)
schubertClassLabels_small = ["\\sigma_{id}", ""]
QH_print_structure_constants_in_basis(F3_schubert_small.self, S_small, schubertClassLabels_small)

# after the output, remove lines which are multiplication by one.