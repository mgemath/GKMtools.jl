using Test
using GKMtools
using Oscar

function is_polynomial_result(x)
  try
    return denominator(x) == 1
  catch
    return true
  end
end

function classification(x)
  return is_polynomial_result(x) ? "polynomial" : "rational"
end

results = String[]

P2 = projective_space(GKM_graph, 2)
beta_P2 = curve_class(P2, Oscar.Edge(1, 2))
P_g0 = ev(1, point_class(P2, 1)) * ev(2, point_class(P2, 2))
gw_default_g0 = gromov_witten(P2, beta_P2, 2, P_g0; show_bar=false)
gw_compact_g0 = gromov_witten(P2, beta_P2, 2, P_g0; show_bar=false, compact_type_only=true)
@test gw_default_g0 == gw_compact_g0
push!(results, "Genus 0 sanity check on P2: compact_type_only leaves the answer unchanged: $gw_compact_g0")

F = gkm_3d_twisted_flag()
beta_F = curve_class(F, "3", "4")
push!(results, "Twisted flag, compact_type_only=true:")
for d in 1:3
  gw = gromov_witten(F, d * beta_F, 0, class_one(); g=1, show_bar=false, compact_type_only=true)
  push!(results, "  genus 1, degree $d: $(classification(gw)) :: $gw")
end

gw_F_d1 = gromov_witten(F, beta_F, 0, class_one(); g=1, show_bar=false, compact_type_only=true)
gw_F_d2 = gromov_witten(F, 2 * beta_F, 0, class_one(); g=1, show_bar=false, compact_type_only=true)
gw_F_d3 = gromov_witten(F, 3 * beta_F, 0, class_one(); g=1, show_bar=false, compact_type_only=true)
@test is_polynomial_result(gw_F_d1)
@test !is_polynomial_result(gw_F_d2)
@test !is_polynomial_result(gw_F_d3)

P1 = projective_space(GKM_graph, 1)
beta_P1 = curve_class(P1, "1", "2")
w = point_class(1, P1)
n = 4
g = 2
d = n + 1 - g
P_psi = prod(i -> ev(i, w), 1:n) * Psi([2 for _ in 1:n])
gw_P1 = gromov_witten(P1, d * beta_P1, n, P_psi; g=g, show_bar=false, compact_type_only=true)
@test is_polynomial_result(gw_P1)
push!(results, "P1 Psi example, compact_type_only=true:")
push!(results, "  genus $g, degree $d, n=$n: $(classification(gw_P1)) :: $gw_P1")

push!(results, "P1 class_one examples, compact_type_only=true:")
for gg in 1:3
  for dd in 1:3
    gw = gromov_witten(P1, dd * beta_P1, 0, class_one(); g=gg, show_bar=false, compact_type_only=true)
    if dd == 1
      @test is_polynomial_result(gw)
    else
      @test !is_polynomial_result(gw)
    end
    push!(results, "  genus $gg, degree $dd: $(classification(gw)) :: $gw")
  end
end

P2_compact = projective_space(GKM_graph, 2)
beta_P2_compact = curve_class(P2_compact, Oscar.Edge(1, 2))
P2_examples = [
  (1, 1, 3, prod(ev(i, point_class(P2_compact, 2)) for i in 1:3)),
  (1, 2, 6, prod(ev(i, point_class(P2_compact, 2)) for i in 1:6)),
  (2, 1, 4, prod(ev(i, point_class(P2_compact, 2)) for i in 1:4)),
  (2, 2, 7, prod(ev(i, point_class(P2_compact, 2)) for i in 1:7)),
]
push!(results, "P2 point-insertion examples, compact_type_only=true:")
for (gg, dd, nn, P_input) in P2_examples
  gw = gromov_witten(P2_compact, dd * beta_P2_compact, nn, P_input; g=gg, show_bar=false, compact_type_only=true)
  if dd == 1
    @test is_polynomial_result(gw)
  else
    @test !is_polynomial_result(gw)
  end
  push!(results, "  genus $gg, degree $dd, n=$nn: $(classification(gw)) :: $gw")
end

G = total_space(vector_bundle_O(1, [-1, -1]))
beta_G = curve_class(G, Oscar.Edge(1, 2))
P_local = Psi(1, 1)
push!(results, "Local Tot(O_P1(-1) \\oplus O_P1(-1)) Psi example, compact_type_only=true:")
for gg in 1:3
  for dd in 1:3
    gw = gromov_witten(G, dd * beta_G, 2, P_local; g=gg, show_bar=false, compact_type_only=true)
    if (gg, dd) in ((1, 1), (1, 2), (1, 3), (2, 1), (3, 1))
      @test is_polynomial_result(gw)
    else
      @test !is_polynomial_result(gw)
    end
    push!(results, "  genus $gg, degree $dd: $(classification(gw)) :: $gw")
  end
end

println(join(results, "\n"))
