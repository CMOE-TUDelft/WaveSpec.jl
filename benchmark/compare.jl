# Compare two result files written by benchmark/benchmarks.jl.
#   julia --project=. benchmark/compare.jl before.json after.json

using BenchmarkTools, Printf

before = minimum(BenchmarkTools.load(ARGS[1])[1])
after  = minimum(BenchmarkTools.load(ARGS[2])[1])

@printf("%-58s %12s %12s %9s %12s %12s\n", "benchmark", "before", "after", "speedup", "mem before", "mem after")
for (key, b) in sort(collect(BenchmarkTools.leaves(before)); by = first)
  a = try after[key] catch; continue end
  @printf("%-58s %12s %12s %8.2fx %12s %12s\n", join(key, " / "),
          BenchmarkTools.prettytime(time(b)), BenchmarkTools.prettytime(time(a)),
          time(b) / time(a), BenchmarkTools.prettymemory(memory(b)), BenchmarkTools.prettymemory(memory(a)))
end
