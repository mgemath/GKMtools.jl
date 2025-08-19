# Case study of the full flag variety F3 for C^3

R = root_system(:A, 2);
F3 = generalized_gkm_flag(R);
BO = get_bruhat_order_of_generalized_flag(R);
S = schubert_classes(F3, BO);
QH_structure_constants(F3);
SB = QH_structure_constants_in_basis(F3, S)
#schubertClassLabels = ["\\sigma_{" * BO.labels[i] * "}" for i in 1:length(BO.labels)];
schubertClassLabels = ["\\sigma_{id}", "\\sigma_{s_1}", "\\sigma_{s_2s_1}", "", "\\sigma_{s_2}", "\\sigma_{s_1s_2}"]
QH_print_structure_constants_in_basis(F3, S, schubertClassLabels)



# after the output, remove lines which are multiplication by one.