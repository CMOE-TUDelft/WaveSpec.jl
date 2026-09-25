using ..AiryWaves
using .PhysicalConstants: g

export AiryRealization, realize

"""
    AiryRealization

A realized Airy-wave field.

Unlike `AiryState`, which stores the spectral and directional
discretization together with a random seed, an `AiryRealization`
stores the fully materialized wave components:

    (ω, kx, ky, amplitude, phase)

and therefore represents one deterministic realization of the sea state.

Two realizations generated from the same `AiryState` will be identical
provided the same random seed is used.

    AiryRealization(components::WaveComponents{T}, k::Vector{T}, h::T)

Build a realization from explicit components, their wavenumber magnitudes `k`
(one per component) and the water depth `h`. The depth-independent factors of
the kinematic kernels are precomputed once per component here, so that the
evaluation kernels only compute the phase and the vertical profile per point:

```math
e_{2kh,i} = e^{-2 k_i h}, \\quad
c_{ϕ,i} =\\frac{g}{\\omega_i (1 + e^{-2 k_i h})}, \\quad
c_{w,i} = \\frac{\\omega_i}{1 - e^{-2 k_i h}}, \\quad
c_{u,i} = c_{w,i} \\frac{k_{x,i}}{k_i}, \\quad
c_{v,i} = c_{w,i} \\frac{k_{y,i}}{k_i}
```
"""
struct AiryRealization{T<:AbstractFloat}
    components::WaveComponents{T}
    k::Vector{T}
    h::T
    e2kh::Vector{T}
    cϕ::Vector{T}
    cu::Vector{T}
    cv::Vector{T}
    cw::Vector{T}

    function AiryRealization(components::WaveComponents{T}, k::Vector{T}, h::T) where {T<:AbstractFloat}
        length(k) == length(components) ||
            throw(DimensionMismatch("k must have one entry per wave component ($(length(components))), got $(length(k))"))
        ω = components.ω
        e2kh = @. exp(-2 * k * h)
        cϕ = @. T(g) / (ω * (1 + e2kh))
        cw = @. ω / (1 - e2kh)
        cu = @. cw * components.kx / k
        cv = @. cw * components.ky / k
        new{T}(components, k, h, e2kh, cϕ, cu, cv, cw)
    end
end

Base.length(realization::AiryRealization) = length(realization.components)

"""
    realize(state::AiryWaves.AiryState)

Build a persistent realization from `state`. The operation only reorganizes
existing Airy-state data and fixes the deterministic random phases associated
with `state.seed`; it does not evaluate the wave field.

Components are stored with the direction index varying fastest: component
`n = (i - 1) * nθ + j` corresponds to frequency `ω[i]` and direction `θ[j]`.
"""
function realize(state::AiryWaves.AiryState)
    nω = state.nω
    nθ = state.nθ

    # get_amplitudes/get_random_phases are (nω × nθ); transpose so that vec()
    # follows the same (i - 1) * nθ + j ordering as ω, k, kx and ky below.
    amplitudes = permutedims(AiryWaves.get_amplitudes(state))
    phases = permutedims(AiryWaves.get_random_phases(state))

    cosθ = cos.(state.θ)
    sinθ = sin.(state.θ)

    ω = repeat(state.ω, inner=nθ)
    k = repeat(state.k, inner=nθ)
    kx = k .* repeat(cosθ, outer=nω)
    ky = k .* repeat(sinθ, outer=nω)

    components = WaveComponents(
        ω, kx, ky, vec(amplitudes), vec(phases)
    )

    return AiryRealization(components, k, state.h)
end
