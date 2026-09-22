using ..AiryWaves

export AiryRealization, realize, evaluate_eta

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
"""
struct AiryRealization{T<:AbstractFloat}
    components::WaveComponents{T}
    k::Vector{T}
    h::T
end

Base.length(realization::AiryRealization) = length(realization.components)

"""
    realize(state::AiryWaves.AiryState)

Build a persistent realization from `state`. The operation only reorganizes
existing Airy-state data and fixes the deterministic random phases associated
with `state.seed`; it does not evaluate the wave field.
"""
function realize(state::AiryWaves.AiryState)
    nω = state.nω
    nθ = state.nθ

    amplitudes = AiryWaves.get_amplitudes(state)
    phases = AiryWaves.get_random_phases(state)

    cosθ = cos.(state.θ)
    sinθ = sin.(state.θ)

    ω = repeat(state.ω, inner=nθ)
    k = repeat(state.k, inner=nθ)
    kx = repeat(state.k, inner=nθ) .* repeat(cosθ, outer=nω)
    ky = repeat(state.k, inner=nθ) .* repeat(sinθ, outer=nω)

    components = WaveComponents(
        ω, kx, ky, vec(amplitudes), vec(phases)
    )

    return AiryRealization(components, k, state.h)
end

"""
    evaluate_eta(realization, x, y, t)

Evaluate free-surface elevation η at a single point and time using a
precomputed Airy realization.

The implementation is intentionally component-wise and allocation-free.
"""
function evaluate_eta(
    realization::AiryRealization{T},
    x::Real,
    y::Real,
    t::Real,
) where {T}

    comps = realization.components

    η = zero(T)

    @inbounds @simd for n in eachindex(comps.ω)

        ψ =
            comps.kx[n] * x +
            comps.ky[n] * y -
            comps.ω[n] * t +
            comps.phase[n]

        η += comps.amplitude[n] * cos(ψ)

    end

    return η
end
