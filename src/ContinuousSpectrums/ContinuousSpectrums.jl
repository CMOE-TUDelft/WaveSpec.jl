module ContinuousSpectrums

using ..Integration
using ..Utils.Truncation

export AbstractSpectrum, JONSWAP, PiersonMoskowitz, TMA, 
       get_density, get_cumulative_energy

# All spectrum models which are defined in this folder/module will inherit from abstract type AbstractSpectrum.
abstract type AbstractSpectrum end

# --- DEFINE AbstractSpectrum's FUNCTIONS INTERFACES ---

"""
    get_density(s::AbstractSpectrum, f::Real)

Returns the spectral density S(f) in [m²s]. 
Must be implemented by all subtypes.
"""
function get_density(s::AbstractSpectrum, f::Real)
    error("Function 'get_density' not implemented for $(typeof(s))")
end


# --- INCLUDE SPECIFIC SPECTRUM MODELS ---

# Include the specific models
include("JONSWAP.jl")
include("TMA.jl")
include("Bretschneider.jl")
include("Donelan.jl")
include("OchiHubble.jl")


# --- SPECTRUM UTILITIES ---

"""
    get_fmax(s::AbstractSpectrum; multiplier=5.0)

Returns a recommended maximum frequency for numerical integration 
based on the peak frequency: f_max = multiplier * f_p.
According to literature, f_max between 3 to 5 times f_p is sufficient to capture more than 99% of the energy.
"""
function get_fmax(s::AbstractSpectrum; multiplier=5.0)
    return multiplier * get_fp(s)
end


"""
    get_cumulative_energy(s::AbstractSpectrum, f, n_points=1000)

Compute the CDF of abstract spectrum.
"""
function get_cumulative_energy(s::AbstractSpectrum, f::Real, n_points::Int=1000)
    fs = range(1e-4, f, length=n_points)
    df = fs[2] - fs[1]
    densities = [get_density(s, fi) for fi in fs]
    return fs, cumsum(densities) .* df
end


end # module