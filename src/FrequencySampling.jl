module FrequencySampling

using ..EnergySpectrum
using ..Integration

export AbstractSampling, UniformSampling, LogSampling, ChebyshevSampling
export generate_grid

abstract type AbstractSampling end

struct UniformSampling   <: AbstractSampling end
struct LogSampling       <: AbstractSampling end
struct ChebyshevSampling <: AbstractSampling end

# --- HELPER: Coordinate Mapping ---
# Maps ξ ∈ [-1, 1] to [a, b]
map_range(ξ, a, b) = 0.5 * (a + b) + 0.5 * (b - a) * ξ

# --- 1. PURE FREQUENCY SAMPLING (No Spectrum Provided) ---
# We map a linear index i in [1, nf] to a frequency f in [fmin, fmax]

function generate_grid(::UniformSampling, fmin, fmax, nf)
    return collect(range(fmin, fmax, length=nf))
end

function generate_grid(::LogSampling, fmin, fmax, nf)
    return exp.(range(log(fmin), log(fmax), length=nf))
end

function generate_grid(::ChebyshevSampling, fmin, fmax, nf)
    # Chebyshev nodes of the first kind in [-1, 1]
    nodes = [cos((2i - 1) * π / (2nf)) for i in nf:-1:1]
    return map_range.(nodes, fmin, fmax)
end

# --- 2. SPECTRAL / ENERGY SAMPLING (Spectrum Provided) ---
# We map a linear index i in [1, nf] to a "cumulative energy fraction" P in [0, 1], 
# then use the Inverse Cumulative Distribution Function (iCDF) of the spectrum to find f.

"""
    generate_grid(strategy, spectrum, fmin, fmax, nf)

Distributes frequencies such that the ENERGY (variance) is spaced 
according to the sampling strategy (e.g., Equal Energy for UniformSampling).
"""
function generate_grid(sampling::AbstractSampling, spec::AbstractSpectrum, fmin, fmax, nf)
    # 1. Compute the Cumulative Distribution Function (CDF)
    # We create a fine-grained lookup table of the integral
    f_lookup = range(fmin, fmax, length=1000)
    df = f_lookup[2] - f_lookup[1]
    
    # Calculate S(f) at all points
    S_vals = [get_density(spec, f) for f in f_lookup]
    
    # Cumulative integration (m0_cumulative)
    m0_cum = cumsum(S_vals) .* df
    m0_total = m0_cum[end]
    
    # Normalize to [0, 1]
    cdf = m0_cum ./ m0_total
    
    # 2. Generate the "Target Energy Probabilities" based on sampling strategy
    target_P = generate_grid(sampling, 0.0, 1.0, nf)
    
    # 3. Inverse Transform: Find frequencies f that match target_P
    # Using linear interpolation to find f for each target probability
    # This is effectively the iCDF
    return [approx_inverse_cdf(target_P[i], cdf, f_lookup) for i in 1:nf]
end

# Simple linear interpolation for the inverse CDF
function approx_inverse_cdf(p_target, cdf, f_lookup)
    p_target <= cdf[1] && return f_lookup[1]
    p_target >= cdf[end] && return f_lookup[end]
    
    # Find the bin where p_target lies
    idx = findfirst(x -> x >= p_target, cdf)
    
    # Linear interpolation
    p0, p1 = cdf[idx-1], cdf[idx]
    f0, f1 = f_lookup[idx-1], f_lookup[idx]
    
    return f0 + (p_target - p0) * (f1 - f0) / (p1 - p0)
end

end # module