# Airy Wave Realizations

WaveSpec distinguishes between:

- `AiryState`
- `AiryRealization`

An `AiryState` defines a stochastic sea state.

```julia
state = AiryState(...)
```

A realization materializes all amplitudes and phases.

```julia
realization = realize(state)
```

The realization can then be evaluated at arbitrary points.

```julia
η = evaluate_η(realization, x, y, t)
```

This representation is useful when repeated evaluations of the same
sea state are required.