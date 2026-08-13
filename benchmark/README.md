# Gromov--Witten benchmarks

Run the benchmarks from the package root:

```sh
julia --project=. benchmark/gromov_witten_fast.jl
julia --project=. benchmark/gromov_witten_nomarks_symbolic.jl
```

Each script performs one untimed warm-up, verifies the invariant is `1`, and
then reports minimum and median wall time, median allocated memory, and median
GC time. No benchmarking dependency is required.

The default is five measured samples. For a quick smoke run or a longer
comparison, set `GKM_BENCH_SAMPLES`:

```sh
GKM_BENCH_SAMPLES=1 julia --project=. benchmark/gromov_witten_fast.jl
GKM_BENCH_SAMPLES=10 julia --project=. benchmark/gromov_witten_nomarks_symbolic.jl
```

Use the same Julia version, thread count, sample count, and machine when
comparing revisions. The minimum time is usually the least noisy signal;
median memory is useful for checking allocation-focused changes.
