# Tests for the twisted quantum product.

V = vector_bundle_O(3, [2])
P3 = baseof(V)
p1 = point_class(P3, 1)
p2 = point_class(P3, 2)
p3 = point_class(P3, 3)
beta = curve_class(P3, "1", "2")

Pred = reduced_virtual_zero_section(V)
P = virtual_zero_section(V)

function tpr(a, b, beta)
    return quantum_product(P3, beta, a, b; useStructureConstants=false, twist_class=Pred, fastMode=false)
end

function tp(a, b, beta)
    return quantum_product(P3, beta, a, b; useStructureConstants=false, twist_class=P, fastMode=false)
end

# Commutativity 1st order:
pr112 = tpr(p1, tpr(p1, p2, beta), beta)
pr121 = tpr(tpr(p1, p2, beta), p1, beta)
# both equal.

p112 = tp(p1, tp(p1, p2, beta), beta)
p121 = tp(tp(p1, p2, beta), p1, beta)
# both equal.

# Associativity 1st order:
p112a = tp(p1, tp(p1, p2, beta), 0*beta) + tp(p1, tp(p1, p2, 0*beta), beta)
p112b = tp(tp(p1, p1, beta), p2, 0* beta) + tp(tp(p1, p1, 0*beta), p2, beta)
# equal

p112a = tpr(p1, tpr(p1, p2, beta), 0*beta) + tpr(p1, tpr(p1, p2, 0*beta), beta)
p112b = tpr(tpr(p1, p1, beta), p2, 0* beta) + tpr(tpr(p1, p1, 0*beta), p2, beta)
# equal.

# Let's try class 2\beta.

pr112a = tp(p1, tp(p1, p2, beta), beta) + tp(p1, tp(p1, p2, 0*beta), 2*beta) + tp(p1, tp(p1, p2, 2*beta), 0*beta)
pr112b = tp(tp(p1, p1, beta), p2, beta) + tp(tp(p1, p1, 0*beta), p2, 2*beta) + tp(tp(p1, p1, 2*beta), p2, 0*beta)
# equal

pr112a = tpr(p1, tpr(p1, p2, beta), beta) + tpr(p1, tpr(p1, p2, 0*beta), 2*beta) + tpr(p1, tpr(p1, p2, 2*beta), 0*beta)
pr112b = tpr(tpr(p1, p1, beta), p2, beta) + tpr(tpr(p1, p1, 0*beta), p2, 2*beta) + tpr(tpr(p1, p1, 2*beta), p2, 0*beta)
# equal

pr123a = tp(p1, tp(p2, p3, beta), beta) + tp(p1, tp(p2, p3, 0*beta), 2*beta) + tp(p1, tp(p2, p3, 2*beta), 0*beta)
pr123b = tp(tp(p1, p2, beta), p3, beta) + tp(tp(p1, p2, 0*beta), p3, 2*beta) + tp(tp(p1, p2, 2*beta), p3, 0*beta)
# NOT EQUAL

pr123a = tpr(p1, tpr(p2, p3, beta), beta) + tpr(p1, tpr(p2, p3, 0*beta), 2*beta) + tpr(p1, tpr(p2, p3, 2*beta), 0*beta)
pr123b = tpr(tpr(p1, p2, beta), p3, beta) + tpr(tpr(p1, p2, 0*beta), p3, 2*beta) + tpr(tpr(p1, p2, 2*beta), p3, 0*beta)
# EQUAL


###
### Example for docs.
###

P3 = projective_space(NormalToricVariety, 3)
L = toric_line_bundle(P3, picard_group(P3)([2]))
l = gkm_line_bundle_of_toric(L)
X = baseof(l)

p1 = point_class(X, 1)
p2 = point_class(X, 2)
beta = curve_class(X, "1", "2")
P = reduced_virtual_zero_section(l)

quantum_product(X, beta, p1, p2; useStructureConstants=false, twist_class=P)