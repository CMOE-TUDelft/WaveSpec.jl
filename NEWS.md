# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.3] - 05/11/2026

### Added
- Realization layer: `WaveComponents`, `AiryRealization`, `realize` and the `evaluate_*` / `evaluate_*!` kernels, type-stable for `Float32` realizations. `AiryRealization(components, k, h)` throws `DimensionMismatch` unless `k` has one entry per component.
- `evaluate_fields` / `evaluate_fields!`: evaluate η, ϕ, u, v and w in a single pass over the wave components (≈4× faster than five separate `evaluate_*` calls).
- Public API for `WaveComponents`, `AiryWaves.get_amplitudes`, and
	`AiryWaves.solve_wavenumber`, with an Airy-realization guide and API reference.
- `benchmark/benchmarks.jl` (BenchmarkTools `SUITE`) and `benchmark/compare.jl` to compare saved benchmark runs. They have their own environment in `benchmark/Project.toml`, so BenchmarkTools is not a dependency of the package.

### Changed
- `generate_sea` now realizes the state once and evaluates the grid point by point with the threaded kernels. Memory use is proportional to the output instead of `Npoints × Ncomponents` (e.g. 841 MiB → 1.8 MiB and 16× faster for a 100×1×5×100 grid with 500 components). Output shapes and keys are unchanged, but the default key order of the returned `NamedTuple` is now `(:η, :ϕ, :u, :v, :w)`, so access fields by name.
- `generate_sea` and the depth profiles in `evaluate_ϕ/u/v/w` use an exact, overflow-free exponential form instead of switching to the deep-water approximation at `kh ≥ 20`. Results for `ϕ`, `u`, `v` and `w` therefore differ slightly from 0.3.2 for components with `kh ≥ 20`, mainly near the bed. `AiryRealization` precomputes the per-component depth-independent factors.

## [0.3.2] - 24/07/2026

### Fixed
- `AiryWaves.change_seed!(state, seed)` threw an error because it read a field that `AiryState` does not have (`state.spec` instead of `state.spectrum`).

## [0.3.1] - 14/04/2026

### Changed
- `WaveSpec.AiryWaves.generate_sea` and `generate_interpolable_sea` can now export the velocity potential field with `vars=[:ϕ]`. Since [PR#36](https://github.com/CMOE-TUDelft/WaveSpec.jl/pull/36).

## [0.3.0] - 29/03/2026

### Changed
- Previous changes were not tracked in the NEWS.md file. Please, see commit history.
