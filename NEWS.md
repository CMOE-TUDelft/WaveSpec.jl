# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `evaluate_fields` / `evaluate_fields!`: evaluate η, ϕ, u, v and w in a single pass over the wave components (≈4× faster than five separate `evaluate_*` calls).
- `benchmark/benchmarks.jl` (BenchmarkTools `SUITE`) and `benchmark/compare.jl` to compare saved benchmark runs.

### Changed
- `generate_sea` now realizes the state once and evaluates the grid point by point with the threaded kernels. Memory use is proportional to the output instead of `Npoints × Ncomponents` (e.g. 841 MiB → 1.8 MiB and 16× faster for a 100×1×5×100 grid with 500 components). Output shapes and keys are unchanged.
- The depth profiles in `evaluate_ϕ/u/v/w` use an exact, overflow-free exponential form instead of switching to the deep-water approximation at `kh ≥ 20`. `AiryRealization` precomputes the per-component depth-independent factors.
- Evaluation kernels are type-stable for `Float32` realizations (previously returned `Union{Float32, Float64}`).
- `AiryRealization(components, k, h)` now throws `DimensionMismatch` unless `k` has one entry per component.
- `WaveSpec.AiryWaves.generate_sea` and `generate_interpolable_sea` can now export the velocity potential field with `vars=[:ϕ]`. Since [PR#36](https://github.com/CMOE-TUDelft/WaveSpec.jl/pull/36).

### Fixed
- `realize` paired amplitudes and phases with the wrong (ω, θ) component when both `nω > 1` and `nθ > 1`, so `evaluate_*` disagreed with `generate_sea` for short-crested seas.

## [0.3.0] - 29/03/2026

### Changed
- Previous changes were not tracked in the NEWS.md file. Please, see commit history.
