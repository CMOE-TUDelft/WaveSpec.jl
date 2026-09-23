using BenchmarkTools
using Random
import WaveSpec.ContinuousSpectrums as CS
import WaveSpec.SpectralSpreading as SS
import WaveSpec.AngularSpreading as AS
import WaveSpec.AiryWaves as AW
import WaveSpec as WS

function create_realization(Ncomponents::Int)
  spec = CS.JONSWAP(1.0, 5.0)
  ds = SS.DiscreteSpectralSpreading(spec,SS.UniformSampling(),1.0,10.0,ceil(Int,Ncomponents/2))
  spread = AS.DiscreteAngularSpreading(0.0) # No angular spreading leads to 2 evaluations
  state = AW.AiryState(ds,spread,50.0)
  realization = WS.realize(state)
  return state, realization
end

function benchmark_evaluation(realization::WS.AiryRealization, Neval::Int; samples::Int = 5, evals::Int = 1)
  rng = Random.Xoshiro(1234)
  
  # Evaluation points
  x = rand(rng, Neval)
  y = rand(rng, Neval)
  t = rand(rng, Neval)
  η = zeros(Float64, Neval)

  # Benchmark
  b = @benchmarkable WS.evaluate_η!($η,$realization,$x,$y,$t,) 
  return run(b; samples=samples, evals=evals)
end

function benchmark_generate_sea(state::WS.AiryState, Neval::Int; samples::Int = 5, evals::Int = 1)
  rng = Random.Xoshiro(1234)
  
  # Evaluation points
  x = rand(rng, Neval)
  y = rand(rng, 1)
  z = Float64[0.0]
  t = rand(rng, 1)

  # Benchmark
  b = @benchmarkable AW.generate_sea($state,$x,$y,$z,$t; vars = [:η]) 
  return run(b; samples=samples, evals=evals)
end

suite = BenchmarkTools.BenchmarkGroup()
for N in (1_000, 10_000, 100_000)
    for Np in (1_000, 10_000, 100_000)
        @info "Benchmark" Ncomponents=N Neval=Np
        state, realization = create_realization(N)
        suite[N, Np]["evaluation"] = benchmark_evaluation(realization, Np)
        suite[N, Np]["generate_sea"] = benchmark_generate_sea(state, Np)
    end
end
display(suite)