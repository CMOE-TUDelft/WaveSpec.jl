# WaveSpec.jl benchmark suite (PkgBenchmark-compatible: defines `SUITE`).
#
# Run and save results:
#   julia --project=. -t auto benchmark/benchmarks.jl results.json
# Compare two saved runs:
#   julia --project=. benchmark/compare.jl before.json after.json

using BenchmarkTools
using Random
import WaveSpec as WS
import WaveSpec.ContinuousSpectrums as CS
import WaveSpec.SpectralSampling as SSa
import WaveSpec.SpectralSpreading as SS
import WaveSpec.AngularSpreading as AS
import WaveSpec.AiryWaves as AW

const NDIR = 10  # directional bins per frequency

"""Short-crested JONSWAP sea with `ncomp ≈ nω × NDIR` components."""
function build_state(ncomp::Int; h = 30.0)
  nω = max(1, ncomp ÷ NDIR)
  ds = SS.DiscreteSpectralSpreading(CS.JONSWAP(2.0, 8.0), SSa.UniformSampling(),
                                    0.05, 0.6, nω + 1; mess = false)
  spread = AS.DiscreteAngularSpreading(:cosinepow, 0.0, 0.5, -π/2, π/2, NDIR + 1)
  return AW.AiryState(ds, spread, h)
end

function points(::Type{T}, n::Int, h) where {T}
  rng = Random.Xoshiro(1234)
  return (T.(100 .* rand(rng, n)), T.(100 .* rand(rng, n)),
          T.(-h .* rand(rng, n)), T.(60 .* rand(rng, n)))
end

const HAS_FUSED = isdefined(WS, :evaluate_fields!)

const SUITE = BenchmarkGroup()

# --- One-off setup cost -------------------------------------------------------
SUITE["setup"] = BenchmarkGroup()
let state = build_state(1000)
  SUITE["setup"]["AiryState_1000"] = @benchmarkable AW.AiryState($(state.spectrum), $(state.spread), 30.0)
  SUITE["setup"]["realize_1000"]   = @benchmarkable WS.realize($state)
end

# --- Point-wise evaluation (hot path) -----------------------------------------
SUITE["evaluate"] = BenchmarkGroup()
for nc in (100, 1000, 4000), np in (1_000, 100_000)
  state = build_state(nc)
  r = WS.realize(state)
  x, y, z, t = points(Float64, np, state.h)
  η, ϕ, u, v, w = (zeros(np) for _ in 1:5)
  g = SUITE["evaluate"]["Nc=$(length(r)) Np=$np"] = BenchmarkGroup()
  g["η"] = @benchmarkable WS.evaluate_η!($η, $r, $x, $y, $t)
  g["ϕ"] = @benchmarkable WS.evaluate_ϕ!($ϕ, $r, $x, $y, $z, $t)
  g["w"] = @benchmarkable WS.evaluate_w!($w, $r, $x, $y, $z, $t)
  # All five fields, one kernel per field (what users must do without a fused API)
  g["all_separate"] = @benchmarkable begin
    WS.evaluate_η!($η, $r, $x, $y, $t)
    WS.evaluate_ϕ!($ϕ, $r, $x, $y, $z, $t)
    WS.evaluate_u!($u, $r, $x, $y, $z, $t)
    WS.evaluate_v!($v, $r, $x, $y, $z, $t)
    WS.evaluate_w!($w, $r, $x, $y, $z, $t)
  end
  if HAS_FUSED
    g["all_fused"] = @benchmarkable WS.evaluate_fields!($η, $ϕ, $u, $v, $w, $r, $x, $y, $z, $t)
  end
end

# --- Single precision ---------------------------------------------------------
SUITE["float32"] = BenchmarkGroup()
let state = build_state(1000), r64 = WS.realize(state), np = 100_000
  c = r64.components
  r = WS.AiryRealization(WS.WaveComponents(Float32.(c.ω), Float32.(c.kx), Float32.(c.ky),
                                           Float32.(c.amplitude), Float32.(c.phase)),
                         Float32.(r64.k), Float32(r64.h))
  x, y, z, t = points(Float32, np, state.h)
  out = zeros(Float32, np)
  SUITE["float32"]["η"] = @benchmarkable WS.evaluate_η!($out, $r, $x, $y, $t)
  SUITE["float32"]["ϕ"] = @benchmarkable WS.evaluate_ϕ!($out, $r, $x, $y, $z, $t)
end

# --- Gridded evaluation through generate_sea ----------------------------------
SUITE["generate_sea"] = BenchmarkGroup()
let state = build_state(500)
  x = collect(range(0.0, 200.0, length = 100)); y = [0.0]
  z = collect(range(-state.h, 0.0, length = 5)); t = collect(range(0.0, 60.0, length = 100))
  SUITE["generate_sea"]["η 100x1x1x100 Nc=500"] = @benchmarkable AW.generate_sea($state, $x, $y, [0.0], $t; vars = [:η])
  SUITE["generate_sea"]["all 100x1x5x100 Nc=500"] = @benchmarkable AW.generate_sea($state, $x, $y, $z, $t)
end

if abspath(PROGRAM_FILE) == @__FILE__
  @info "Running WaveSpec benchmarks" threads = Threads.nthreads() fused = HAS_FUSED
  results = run(SUITE; verbose = true, samples = 5, evals = 1, seconds = 10)
  display(results)
  isempty(ARGS) || BenchmarkTools.save(ARGS[1], results)
end
