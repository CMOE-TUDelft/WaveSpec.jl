"""
    WaveComponents

Structure-of-arrays representation of discrete linear wave components.
"""
struct WaveComponents{T<:AbstractFloat}
    ω::Vector{T}
    kx::Vector{T}
    ky::Vector{T}
    amplitude::Vector{T}
    phase::Vector{T}

    function WaveComponents(
        ω::Vector{T}, kx::Vector{T}, ky::Vector{T},
        amplitude::Vector{T}, phase::Vector{T},
    ) where {T<:AbstractFloat}
        n = length(ω)
        length(kx) == n || throw(DimensionMismatch("kx must have the same length as ω"))
        length(ky) == n || throw(DimensionMismatch("ky must have the same length as ω"))
        length(amplitude) == n || throw(DimensionMismatch("amplitude must have the same length as ω"))
        length(phase) == n || throw(DimensionMismatch("phase must have the same length as ω"))
        new{T}(ω, kx, ky, amplitude, phase)
    end
end

Base.length(components::WaveComponents) = length(components.ω)