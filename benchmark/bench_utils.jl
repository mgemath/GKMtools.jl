module BenchUtils

using Printf

export run_benchmark

"""Run `f` after one untimed warm-up and print stable, dependency-free metrics."""
function run_benchmark(name::AbstractString, f; expected, samples::Int=5)
  samples > 0 || throw(ArgumentError("samples must be positive"))

  warmup = f()
  warmup == expected || error("$name returned $warmup; expected $expected")

  times = Float64[]
  bytes = Int[]
  gc_times = Float64[]

  for _ in 1:samples
    GC.gc()
    measurement = @timed f()
    measurement.value == expected ||
      error("$name changed result during benchmarking: $(measurement.value)")
    push!(times, measurement.time)
    push!(bytes, measurement.bytes)
    push!(gc_times, measurement.gctime)
  end

  sorted_times = sort(times)
  median_time = sorted_times[cld(length(sorted_times), 2)]
  median_bytes = sort(bytes)[cld(length(bytes), 2)]

  println(name)
  @printf("  samples:       %d\n", samples)
  @printf("  minimum time:  %.6f s\n", minimum(times))
  @printf("  median time:   %.6f s\n", median_time)
  @printf("  median memory: %.3f MiB\n", median_bytes / 2.0^20)
  @printf("  median GC:     %.6f s\n", sort(gc_times)[cld(length(gc_times), 2)])
  println("  result:        ", warmup)

  return (; minimum_time=minimum(times), median_time, median_bytes,
          result=warmup)
end

end
