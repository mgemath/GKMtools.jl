using GKMtools
using Oscar

include("bench_utils.jl")
using .BenchUtils

# Symbolic no-marks invariant. This exercises edge-integral evaluation and
# fraction arithmetic rather than the numeric specialization used by fast mode.
P2 = projective_space(GKMGraph, 2)
beta = curve_class(P2, Oscar.Edge(1, 2))
point = point_class(P2, 1)
insertions = [point for _ in 1:5]

samples = parse(Int, get(ENV, "GKM_BENCH_SAMPLES", "5"))

run_benchmark(
  "gromov_witten_nomarks: P2 degree 2, five point insertions, symbolic",
  () -> gromov_witten_nomarks(
    P2, 2 * beta, insertions;
    show_bar=false, fast_mode=false,
  );
  expected=QQ(1), samples,
)
#=
gromov_witten_nomarks: P2 degree 2, five point insertions, symbolic
  samples:       5
  minimum time:  0.249119 s
  median time:   0.257563 s
  median memory: 2.169 MiB
  median GC:     0.000000 s
  result:        1
(minimum_time = 0.249118656, median_time = 0.257562771, median_bytes = 2274400, result = 1)
=#