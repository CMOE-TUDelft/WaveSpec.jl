# WaveSpec.jl performance and code audit: `wave_components` branch

**Scope:** `wave_components` at `06f7a53`, focused on the new realization/evaluation layer (`WaveComponents`, `AiryRealization`, `evaluate_*`) and the older `generate_sea` path it is meant to replace.
**Fixes:** local branch `perf/wave-components-audit` (uncommitted).
**Machine:** Apple Silicon, 14 cores (10 performance), Julia 1.11.1, 8 threads unless stated.

## 1. Summary

| # | Finding | Severity | Status |
|---|---|---|---|
| C1 | `realize` pairs amplitudes and phases with the wrong (ω, θ) component for short-crested seas | **Critical** (wrong results) | Fixed and tested |
| C2 | `generate_sea` materialises 6-D `Np × Nc` tensors: memory is O(points × components) | **High** | Fixed: 6.5–16.6× faster, −99.7% memory |
| C3 | Out-of-bounds reads under `@inbounds` when `length(k) ≠ length(components)` | **High** (memory safety) | Fixed: constructor validates |
| P1 | Computing all fields costs 5 separate passes (5× phase, trig, hyperbolics) | Medium | Fixed: fused `evaluate_fields` gives ~4.5× |
| P2 | Float32 realizations return `Union{Float32, Float64}` and compute in Float64 | Medium | Fixed: inferred `Float32` |
| P3 | Setup (`AiryState`, `realize`) re-samples angles and rebuilds grids in every getter; abstract fields | Medium for ensembles | Recommendation |
| A1 | `kh < 20` branch: discontinuous deep-water switch (up to 2× relative error at the bed) and a branch in the hot loop | Low | Fixed: exact, overflow-free form |
| P4 | `@simd` has no effect: Base `sin/cos/exp` do not vectorise; trig dominates (~12 ns per component per point) | Info | Recommendation (§5) |
| P5 | `@inbounds Threads.@threads …` does not propagate into the threaded closure | Low | Fixed: shared driver |
| E1 | Julia 1.11.1: multi-threaded runs calling `GC.gc()` (BenchmarkTools) can deadlock | Environment | Workaround documented |

Threading is **not** a bottleneck: point-parallel evaluation scales near-linearly up to the 10 performance cores (§3.3).

## 2. Findings in detail

### C1: `realize` mis-pairs components (critical)

`get_amplitudes(state)` and `get_random_phases(state)` are `(nω × nθ)` matrices, so `vec` runs **ω fastest**. `realize` builds `ω`, `k`, `kx`, `ky` with `repeat(…, inner=nθ)`, which runs **θ fastest**. When `nω > 1` and `nθ > 1`, each wave gets another component's amplitude and phase. Reproduction with a JONSWAP spectrum, 5 frequencies × 4 directions and h = 30 m:

```
point (3, -2, -4, 1.5):  generate_sea η = 0.7728   evaluate_η = 0.2711
                         generate_sea u = 0.5068   evaluate_u = 0.3355
```

The existing tests missed this because every realization test used either one frequency (`RegularWave`) or one direction, and `AiryRealizationTests` asserted the inconsistent layout. **Fix:** `permutedims` before `vec`. New tests check the pairing element by element, and check all five fields against an independent textbook double sum over (ω, θ) for a short-crested sea.

### C2: `generate_sea` memory explosion (high)

The phase tensor `ψ` has shape `(nx, ny, nz, nt, nω, nθ)`, and each field creates further temporaries of the same size (`A .* cos.(ψ)`, …). Memory therefore scales as points × components. The benchmark grid (5·10⁴ points, 500 components) already allocates 841 MiB. A modest CFD inflow case (100×100×10 points × 1000 time steps × 1000 components) would need about 800 GB per temporary.

**Fix:** `generate_sea` realizes the state once and evaluates the grid point by point with the threaded, allocation-free kernels (`η` on `(nx, ny, 1, nt)`, fused kinematics on `(nx, ny, nz, nt)`). Output shapes, keys and key order are unchanged. The method is declared in `AiryWaves` (`function generate_sea end`) and implemented in `AiryEvaluation.jl`, because it needs `realize`, which is defined after `AiryWaves`.

### C3: unchecked `k` length (high, memory safety)

`AiryRealization(components, k, h)` accepted any `k`, and the kernels index `k[i]` under `@inbounds`. The branch's own tests built realizations with 3 components and `k = [0.0]`, so any `evaluate_ϕ/u/v/w` call on them reads out of bounds. The constructor now throws `DimensionMismatch`, and the in-place drivers use `eachindex(out, x, y, …)`. That call checks that all axes match, which makes `@inbounds` sound.

### P1: redundant work across fields (medium)

A user who needs η, ϕ, u, v and w (the usual case for wave-maker or inflow boundary conditions) previously paid for five phase evaluations, five trig calls and six `cosh`/`sinh` per component per point. The new `evaluate_fields` and `evaluate_fields!` compute `sincos(ψ)` and the vertical profile once, then accumulate all five sums. Their cost is 1.06–1.23× that of a single `evaluate_ϕ`.

### P2: Float32 type instability (medium)

The literal `20.0` and the Float64 constant `g` promoted every Float32 kernel to Float64 on some paths: `Base.return_types(evaluate_ϕ, (AiryRealization{Float32}, Float32, …))` gave `Union{Float32, Float64}`. Now `g` is converted once in the constructor (`T(g)`), and the accumulator type is `promote_type(T, typeof.(coords)...)`. `@inferred` tests cover `evaluate_ϕ`, `evaluate_w` and `evaluate_fields`. Scalar Float32 is **not** faster on CPU, since the trig costs about the same. Its value is memory and future GPU use.

### P3: setup path (medium for ensembles)

`realize` takes ~7 ms and `AiryState` ~2 ms for 1000 components (10 directions):
- `get_central_angles`, `get_bandwidths` and `get_weights` each call `get_angles`, which re-samples and sorts `nθ` random angles every time.
- `get_frequencies` regenerates the grid. With `Energy` sampling this means repeated inverse-CDF integrations.
- `AiryState.spectrum::DiscreteSpectralSpreading` and `AiryState.spread::DiscreteAngularSpreading` are UnionAll (abstract) fields, and `DiscreteSpectralSpreading` has `fmin/fmax/norm_factor::Real` and `domain::SamplingDomain`. As a result `AiryWaves.get_amplitudes(state)` infers as `Any`.

This is harmless for one sea state but dominates for Monte Carlo ensembles: 1000 seeds cost about 9 s of setup before any evaluation. See §5 for the fix.

### A1: deep-water switch (low, accuracy)

The kernels used `cosh(k(z+h))/cosh(kh)` for `kh < 20` and `exp(kz)` otherwise. The switch is discontinuous: at the bed the exact ratio is `1/cosh(kh)` while the approximation gives `e^{-kh}`, a factor of about 2. Absolute values are tiny, but the branch sat inside the hot loop. The replacement is exact for all `kh` and cannot overflow:

```math
\frac{\cosh k(z+h)}{\cosh kh} = \frac{e^{kz} + e^{-k(z+2h)}}{1 + e^{-2kh}},\qquad
\frac{\cosh k(z+h)}{\sinh kh} = \frac{e^{kz} + e^{-k(z+2h)}}{1 - e^{-2kh}},\qquad
\frac{\sinh k(z+h)}{\sinh kh} = \frac{e^{kz} - e^{-k(z+2h)}}{1 - e^{-2kh}}
```

`e^{-2kh}` and the denominators are precomputed per component in `AiryRealization`, and `e^{-k(z+2h)} = e^{-2kh}/e^{kz}`. That leaves **one** `exp` per component per point, 19% cheaper than two, with a branchless guard for the case where both terms underflow. A test with `kh > 700` (where `cosh(kh)` overflows) checks that results are finite and match the deep-water limit near the surface.

## 3. Benchmarks

The suite lives in `benchmark/benchmarks.jl`. It defines a PkgBenchmark-compatible `SUITE` with these groups:
- `setup`
- `evaluate`: Nc ∈ {100, 1000, 4000} × Np ∈ {10³, 10⁵}, for η, ϕ, w, all fields with separate kernels, and all fields fused
- `float32`
- `generate_sea`

`benchmark/compare.jl` prints before/after tables. The baseline was measured on a clean worktree at `06f7a53`. Times are minima over 5 samples.

### 3.1 Before vs after

| Benchmark | Before | After | Speed-up | Memory before → after |
|---|---|---|---|---|
| `generate_sea` η, 100×1×1×100, Nc=500 | 97.2 ms | 14.9 ms | **6.5×** | 76.5 MiB → 244 KiB (−99.7%) |
| `generate_sea` all fields, 100×1×5×100, Nc=500 | 1.012 s | 61.1 ms | **16.6×** | 841 MiB → 1.8 MiB (−99.8%) |
| all 5 fields, Nc=1000, Np=10⁵ (separate → fused) | 1.003 s | 207 ms | **4.9×** | 24 KiB → 5 KiB |
| all 5 fields, Nc=4000, Np=10⁵ (separate → fused) | 4.181 s | 1.112 s | **3.8×** | 24 KiB → 5 KiB |
| `evaluate_ϕ!`, Nc=1000, Np=10⁵ | 220 ms | 189 ms | 1.17× | unchanged |
| `evaluate_η!`, Nc=1000, Np=10⁵ | 156 ms | 162 ms | 0.96× (noise) | unchanged |
| `realize`, Nc=1000 | 6.5 ms | 7.2 ms | 0.91× | 187 → 226 KiB (precomputed coefficients) |

Fused kernel against separate calls, across sizes:

| Size | 5 separate calls (before) | `evaluate_fields!` | Speed-up | Fused / single `ϕ` |
|---|---|---|---|---|
| Nc=100, Np=10³ | 1.04 ms | 211 μs | 4.9× | 1.12 |
| Nc=100, Np=10⁵ | 94.6 ms | 19.6 ms | 4.8× | 1.08 |
| Nc=1000, Np=10³ | 10.1 ms | 2.16 ms | 4.7× | 1.06 |
| Nc=1000, Np=10⁵ | 1.00 s | 207 ms | 4.9× | 1.10 |
| Nc=4000, Np=10³ | 41.4 ms | 10.3 ms | 4.0× | 1.21 |
| Nc=4000, Np=10⁵ | 4.18 s | 1.11 s | 3.8× | 1.23 |

Single-field kernels are unchanged within run-to-run variation (±7%). They are bound by `sin/cos`, and the new form mainly buys exactness, robustness and type stability. The in-place drivers allocate about 5 KiB per call, which is the task-spawn overhead of `Threads.@threads`. The scalar kernels allocate nothing, and a test asserts this.

### 3.2 Numerical accuracy

- All 729 package tests pass, up from 352; the new ones cover C1–C3, P1, P2 and A1.
- `evaluate_*`, `evaluate_fields` and `generate_sea` agree with an independent textbook (ω, θ) double sum to `atol = 1e-12` on a short-crested sea, at points down to the bed.
- Float32 results match Float64 to `rtol = 1e-4`.

### 3.3 Thread scaling (Nc = 1000, Np = 2·10⁴)

| Threads | η before | η after | ϕ before | ϕ after | fused after |
|---|---|---|---|---|---|
| 1 | 237 ms | 246 ms | 335 ms | 289 ms | 314 ms |
| 2 | 120 ms | 124 ms | 170 ms | 146 ms | 157 ms |
| 4 | 61 ms | 62 ms | 86 ms | 74 ms | 79 ms |
| 8 | 31 ms | 32 ms | 46 ms | 40 ms | 44 ms |
| 10 | 29 ms | 28 ms | 41 ms | 35 ms | 38 ms |

Parallel efficiency is about 95% at 8 threads. Even Np = 100 scales (1.0 ms → 0.1 ms on 10 threads), so no serial cutoff is needed at these sizes.

## 4. Changes made (branch `perf/wave-components-audit`)

- `src/AiryRealization.jl`:
  - fixed the component ordering in `realize`
  - `AiryRealization` validates `k` and precomputes `e2kh`, `cϕ`, `cu`, `cv`, `cw`
- `src/AiryEvaluation.jl`:
  - exact single-`exp` vertical profiles
  - type-generic accumulators
  - shared `_pointwise!` threaded driver with axis checking
  - new `evaluate_fields` / `evaluate_fields!`
  - the low-memory `generate_sea` method
- `src/AiryWaves.jl`: the tensor implementation of `generate_sea` is replaced by a documented function stub.
- Tests:
  - pairing tests
  - textbook-reference tests for a short-crested sea
  - deep-water (`kh > 700`) test
  - `evaluate_fields!` and Float32 inference tests
  - zero-allocation test
  - existing hand-built realizations now pass a `k` of matching length
- `benchmark/benchmarks.jl`, `benchmark/compare.jl`, and a `NEWS.md` entry.

Not verified locally: **Julia 1.10**, which CI tests. The code avoids newer syntax, but this should be confirmed in CI.

## 5. Recommendations

**Performance (ordered by expected gain):**

1. **Phase recurrence for structured grids (4–5×, measured).** For a uniform time step, `e^{iψ(t+Δt)} = e^{iψ(t)} e^{-iωΔt}`, so each step needs one complex multiply per component instead of `sincos`. A scratch prototype for a η time series (Nc = 1000, 2·10⁴ steps, one thread), re-seeded with an exact `sincos` every 256 steps, ran **4.6× faster** (164 → 36 ms) with a maximum error of 2·10⁻¹³ m. The same idea applies along uniform `x`/`y` grid lines. It would directly benefit `Signal.generate_signal`, `generate_sea` on regular grids, and wave-maker time series. Expose it as `evaluate_*_timeseries!` or as a `generate_sea` fast path when the inputs are `AbstractRange`s.
2. **Vectorise across components.** The loop is scalar because Base trig does not SIMD-vectorise. A branch-free polynomial `sincos` for already-reduced arguments, over SIMD.jl lanes or via a SLEEF-style library, could give 2–4× on the remaining exact path. Avoid LoopVectorization.jl: it is not maintained for Julia ≥ 1.11.
3. **GPU extension.** The struct-of-arrays `WaveComponents` layout and the pure scalar kernels map directly to one GPU thread per point looping over components. A `KernelAbstractions.jl` package extension (weak dependency) would keep the core CPU-only. This is where Float32 pays off.
4. **Ensembles.** Split the deterministic part (ω, k, θ, amplitudes, precomputed coefficients) from the random phases, so that `realize(state, seed)` for a new seed only redraws `Nc` phases. Compute the angles once in `DiscreteAngularSpreading`, or cache them, instead of in every getter.

**Type stability and structure:**

5. Parametrise `AiryState{S<:DiscreteSpectralSpreading, D<:DiscreteAngularSpreading}` and give `DiscreteSpectralSpreading` concrete `fmin/fmax/norm_factor` and a `domain` type parameter. The setup path will then infer concretely.
6. Rename `change_seed!` to `with_seed` (or make it actually mutate). It returns a new immutable state, and callers who ignore the return value silently keep the old seed.
7. Move `WaveComponents`, `AiryRealization` and the evaluation kernels into a submodule (for example `AiryWaves.Realization`) included **before** `generate_sea`'s consumers. That removes the cross-module method definition used here and lets `@autodocs` pick up the new API, which is currently missing from the docs. Decide on exports consistently: the `evaluate_*!` functions are not exported, while their scalar counterparts are.
8. Keep `generate_sea` as a thin compatibility wrapper over the realization API, and consider returning a concretely typed `NamedTuple` (for example always all 5 fields, or a `vars::Val` argument) for type-stable callers.

**Housekeeping:**

9. Move `BenchmarkTools` and `Test` out of `[deps]`, into `[extras]`/test targets and a `benchmark/Project.toml`. Retire `benchmark/AiryEvaluationBenchmarks.jl`, which is superseded, uses only one direction, and has a comment claiming "2 evaluations".
10. **Environment (E1):** on Julia 1.11.1, BenchmarkTools runs with `-t 8` hung twice with all worker threads asleep and the main thread in `jl_gc_wait_for_the_world` under an explicit `GC.gc()`. The baseline code hung the same way, so this is a runtime issue, not WaveSpec's. `JULIA_THREAD_SLEEP_THRESHOLD=infinite` avoided it. Upgrade to the latest 1.11.x patch release and re-check.
