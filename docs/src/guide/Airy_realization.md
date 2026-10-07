# Airy Wave Realizations

An `AiryState` stores the frequency and direction discretization together with
a seed. `realize` fixes the phases and materializes the wave components for
repeated pointwise evaluation.

```julia
using WaveSpec
import WaveSpec.AiryWaves as AW

spec = WaveSpec.ContinuousSpectrums.RegularWave(0.4, 6.0)
discrete = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
	spec; mess = false,
)
spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
state = AW.AiryState(discrete, spread, 12.0)
realization = realize(state)

fields = evaluate_fields(realization, 1.2, 0.0, -2.0, 1.0)
η = evaluate_η(realization, 1.2, 0.0, 1.0)
```

`AiryRealization` stores a `WaveComponents` object with flat `ω`, `kx`, `ky`,
amplitude, and phase vectors, plus one wavenumber magnitude per component.
Components produced by `realize` are ordered with direction varying fastest:
`n = (i - 1) * nθ + j` maps frequency `i` and direction `j` to component `n`.

```julia
components = WaveComponents(
	[1.0], [0.1], [0.0], [0.2], [0.0],
)
explicit = AiryRealization(components, [0.11], 12.0)
evaluate_fields(explicit, 1.2, 0.0, -2.0, 1.0)
```

For `evaluate_ϕ`, `evaluate_u`, `evaluate_v`, and `evaluate_w`, `z` is measured
upward from the still-water surface at `z = 0`; `evaluate_η` is independent of
`z`. Scalar evaluators return values at one point. `evaluate_fields` computes
all five fields in one pass. The in-place, threaded `evaluate_*!` methods are
available through qualified names such as `WaveSpec.evaluate_fields!`.

For the discrete amplitude and phase matrices of an `AiryState`, use
`AW.get_amplitudes(state)` and `AW.get_random_phases(state)`; both have shape
`(nω, nθ)`. `AW.solve_wavenumber(ω, h)` solves the finite-depth dispersion
relation used to construct explicit realizations.
