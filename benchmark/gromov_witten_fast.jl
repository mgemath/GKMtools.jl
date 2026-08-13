using GKMtools
using Oscar

include("bench_utils.jl")
using .BenchUtils

# Enumeration-heavy marked invariant. The warm-up is deliberately excluded so
# compilation and package loading do not obscure changes in the localization
# graph traversal.
P2 = projective_space(GKMGraph, 2)
beta = curve_class(P2, Oscar.Edge(1, 2))
point = point_class(P2, 1)
n_marks = 5
insertion = prod(i -> ev(i, point), 1:n_marks)

samples = parse(Int, get(ENV, "GKM_BENCH_SAMPLES", "5"))

run_benchmark(
  "gromov_witten: P2 degree 2, five marked points, fast mode",
  () -> gromov_witten(
    P2, 2 * beta, n_marks, insertion;
    show_bar=false, fast_mode=true,
  );
  expected=QQ(1), samples,
)
#=
gromov_witten: P2 degree 2, five marked points, fast mode
  samples:       5
  minimum time:  0.250270 s
  median time:   0.261225 s
  median memory: 5.302 MiB
  median GC:     0.000000 s
  result:        1
(minimum_time = 0.250269967, median_time = 0.261225468, median_bytes = 5559784, result = 1)

=#