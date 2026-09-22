using ..AiryWaves

export AiryRealization, realize

"""
    AiryRealization

Persistent realization of the discrete Airy-wave components.

This is the first step toward matrix-free wave-field evaluation: the
frequency-direction pair `(i,j)` is flattened into one component index.
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
