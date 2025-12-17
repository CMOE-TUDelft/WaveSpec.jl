module SpectralSpreading

using ..FrequencySampling
using ..EnergySpectrum
using ..AngularSpreading

export Spectrum, discretize_spectrum

struct Spectrum
    frequencies::Vector{Float64}
    densities::Vector{Float64}     # S(f)
    angles::Vector{Vector{Float64}} # Vector of vectors if f-dependent
    weights::Vector{Vector{Float64}}
end

"""
    discretize_spectrum(spec_shape, sampling_strat, f_min, f_max, n_f)
    
The master function that coordinates the Grid and the Shape.
"""
function discretize_spectrum(shape::AbstractSpectrum, strat::AbstractSampling, f_min, f_max, n_f)
    # 1. Logic for Grid
    freqs = if strat isa EqualEnergySampling
        # Special logic calling EnergySpectrum.get_cumulative_energy
    else
        FrequencySampling.generate_grid(strat, f_min, f_max, n_f)
    end
    
    # 2. Map densities
    densities = [get_density(shape, f) for f in freqs]
    
    return freqs, densities
end

end # module