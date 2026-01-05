module FrequencySampling

using ..ContinuousSpectrums
using ..Integration

export AbstractSampling, UniformSampling, LogSampling, ChebyshevSampling
export SamplingDomain, FrequencyDomain, EnergyDomain, Frequency, Energy
export generate_grid

# ---- Define sampling strategies markers ----

abstract type AbstractSampling end

struct UniformSampling   <: AbstractSampling end
struct LogSampling       <: AbstractSampling end
struct ChebyshevSampling <: AbstractSampling end

# --------------------------------------------

# ---- Define sampling domain markers ---- 

abstract type SamplingDomain end

struct FrequencyDomain <: SamplingDomain end
struct EnergyDomain    <: SamplingDomain end

# Constants for easy use: FrequencySampling.Frequency
const Frequency = FrequencyDomain()
const Energy    = EnergyDomain()

# ----------------------------------------

# --- HELPER: Coordinate Mapping ---
# Maps ξ ∈ [-1, 1] to [a, b]
map_range(ξ, a, b) = 0.5 * (a + b) + 0.5 * (b - a) * ξ

# --- BASE SAMPLING ( EQUIVALENT TO SAMPLE ON FREQUENCY DOMAIN) ---
# We map a linear index i in [1, nf] to a frequency f in [fmin, fmax]

function generate_grid(::UniformSampling, fmin, fmax, nf)
    return collect(range(fmin, fmax, length=nf))
end

function generate_grid(::LogSampling, fmin, fmax, nf)
    # If fmin is 0, we must nudge it to avoid -Inf
    _fmin = fmin <= 0.0 ? 1e-6 : fmin
    return exp.(range(log(_fmin), log(fmax), length=nf))
end

function generate_grid(::ChebyshevSampling, fmin, fmax, nf)
    # Chebyshev nodes of the first kind in [-1, 1]
    nodes = [cos((2i - 1) * π / (2nf)) for i in nf:-1:1]
    return map_range.(nodes, fmin, fmax)
end

# Now we implement a common signature (AbstractSampling, SamplingDomain, AbstractSpectrum...) method for 
# both frequency and energy sampling which will depend on the Domain marker argument. 

#   generate_grid(sampling::AbstractSampling, ::SamplingDomain, spec::AbstractSpectrum, fmin, fmax, nf)

# For frequency domain, the sampling is equivalent to the base sampling methods defined above.

# --- FREQUENCY (BASE) SAMPLING ---
"""
    generate_grid(strategy, domain, spectrum, fmin, fmax, nf)

Distributes frequencies values according to the sampling strategy (UniformSampling, LogSampling...).
"""
function generate_grid(sampling::AbstractSampling, ::FrequencyDomain, spec::AbstractSpectrum, fmin, fmax, nf)
    return generate_grid(sampling, fmin, fmax, nf) 
end

# --- ENERGY SPECTRAL SAMPLING -------
# For energy sampling, we map a linear index i in [1, nf] to a "cumulative energy fraction" P in [0, 1], 
# then use the Inverse Cumulative Distribution Function (iCDF) of the spectrum to find f.
"""
    generate_grid(strategy, domain, spectrum, fmin, fmax, nf)

Distributes frequencies such that the ENERGY (variance) is spaced 
according to the sampling strategy (e.g., Equal Energy for UniformSampling).
"""
function generate_grid(sampling::AbstractSampling, ::EnergyDomain, spec::AbstractSpectrum, fmin, fmax, nf)
    # 1. Compute energy CDF lookup
    freqs_fine = range(fmin, fmax, length=1000)
    central_freqs_fine = (freqs_fine[1:end-1] + freqs_fine[2:end]) ./ 2.0
    Δf_fine = diff(freqs_fine)
    S_fine = [get_density(spec, f) for f in central_freqs_fine]
    cdf = vcat(0.0, cumsum(S_fine .* Δf_fine))
    cdf ./= cdf[end] # Normalize so that total energy = 1
    
    # 2. Get target normalised energy values using the sampling strategy on range (0, 1)
    # e.g., if sampling is Uniform, Ei is linspace(0, 1)
    Ei = generate_grid(sampling, 1.0e-5, 1.0, nf)
    
    # 3. Apply the optimized iCDF
    return approx_inverse_cdf(cdf, freqs_fine, Ei, nf)
end

"""
    approx_inverse_cdf(p_target, cdf, f_lookup)

Uses binary search and linear interpolation to find the frequency `f` 
corresponding to a cumulative probability `p_target`.
"""
function approx_inverse_cdf(cdf, f_fine, Ei, nf)
    f_target = zeros(nf)

    for (i, E_target) in enumerate(Ei)
        # Handle boundaries explicitly to avoid idx=1 or idx > length
        if E_target <= cdf[1]
            f_target[i] = f_fine[1]
            continue
        elseif E_target >= cdf[end]
            f_target[i] = f_fine[end]
            continue
        end

        # Binary Search
        idx = searchsortedfirst(cdf, E_target)
        
        # Linear Interpolation
        p0, p1 = cdf[idx-1], cdf[idx]
        f0, f1 = f_fine[idx-1], f_fine[idx]
    
        f_target[i] = f0 + (E_target - p0) * (f1 - f0) / (p1 - p0)
    end

    return f_target
end

end # module