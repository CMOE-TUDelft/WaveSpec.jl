# Getting Started

This guide covers the basic workflow of using WaveSpec.jl to generate stochastic sea states.

## Basic Workflow

The typical workflow involves four steps:
1.  Defining a **Continuous Spectrum**.
2.  **Discretizing** the spectrum into frequency components.
3.  Defining **Angular Spreading**.
4.  Generating the **Airy State**.

### 1. Define a Continuous Spectrum

Choose a spectral model (e.g., `JONSWAP`, `Bretschneider`, `TMA`).

```julia
using WaveSpec

Hs = 2.5  # Significant wave height [m]
Tp = 8.0  # Peak period [s]
spec = JONSWAP(Hs, Tp)
```

### 2. Discretize the Spectrum

Convert the continuous spectrum into discrete frequency bands. You can choose different sampling strategies like `LogSampling` or `UniformSampling`.

```julia
# DiscreteSpectralSpreading(shape, strategy, fmin, fmax, nf; domain=Frequency)
discrete_spec = DiscreteSpectralSpreading(
    spec,
    LogSampling(),         # Sampling strategy
    0.05,                  # fmin [Hz]
    1.0,                   # fmax [Hz]
    50;                    # Number of frequencies
    domain = Frequency     # Sampling domain
)
```

### 3. Define Angular Spreading

Define how the energy is distributed across directions.

```julia
# Cosine-Power distribution
angle_dist = CosinePowerDistribution(0.0, 2.0)

# Discretize into 36 angles between -π/2 and π/2
discrete_angle = DiscreteAngularSpreading(angle_dist, -π/2, π/2, 36)
```

### 4. Given the Sea State

Combine the spectral and angular information with the water depth.

```julia
depth = 50.0 # Water depth [m]
sea_state = AiryState(discrete_spec, discrete_angle, depth)
```

You can now access the wave components:
- `sea_state.ω`: Angular frequencies
- `sea_state.k`: Wavenumbers
- `sea_state.θ`: Directions
