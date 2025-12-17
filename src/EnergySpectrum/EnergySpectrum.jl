module EnergySpectrum

using ..Integration

export AbstractSpectrum, JONSWAP, PiersonMoskowitz, TMA,
       get_density, get_cumulative_energy, estimate_γ

abstract type AbstractSpectrum end

# Include the specific models
include("JONSWAP.jl")
include("TMA.jl")

"""
    get_cumulative_energy(s::AbstractSpectrum, ω_max, n_points=1000)
General utility to find the CDF of any spectrum for Equal Energy sampling.
"""
function get_cumulative_energy(s::AbstractSpectrum, ω_max::Real, n_points::Int=1000)
    ωs = range(1e-4, ω_max, length=n_points)
    dω = ωs[2] - ωs[1]
    densities = [get_density(s, ω) for ω in ωs]
    
    # Returns the frequency axis and the cumulative sum (m0 development)
    return ωs, cumsum(densities) .* dω
end

end # module