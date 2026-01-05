module SpectralSpreading

using ..FrequencySampling
using ..ContinuousSpectrums
using ..Integration

export DiscreteSpectrum
export get_frequency, get_frequencies,
       get_density, get_densities,
       get_amplitude, get_amplitudes 
export get_moments, get_Hs, discrete_cdf

"""
    DiscreteSpectrum
The discrete realization of a frequency spectrum. 
Represents the energy distribution across a finite set of frequency bins.
"""
struct DiscreteSpectrum
    spectrum::AbstractSpectrum  # The underlying continuous spectrum (JONSWAP, etc.)
    sampling::AbstractSampling  # The sampling strategy used
    domain::SamplingDomain      # Domain over which the spectral sampling is done (Frequency/Energy)
    fmin::Real                  # Lower frequency [Hz]
    fmax::Real                  # Upper frequency [Hz]
    nf::Int                     # Number of frequency samples
    nbands::Int                 # Number of frequency bands ( = nf -1 )
    norm_factor::Real           # Normalization factor applied
end

# --- CONSTRUCTOR ---

"""
    DiscreteSpectrum(shape, strategy, fmin, fmax, nf)

Continuous abstract spectrum discretised according to the selected abstract sampling strategy, between fmin and fmax with nf bins.
"""
function DiscreteSpectrum(shape::AbstractSpectrum, strat::AbstractSampling, fmin::Real, fmax::Real, nf::Int; domain::SamplingDomain = Frequency, mess::Bool=true)
    # 1. Generate the nf edges
    freqs = FrequencySampling.generate_grid(strat, domain, shape, fmin, fmax, nf)
    
    # 2. Compute bandwidths and centers for the nf-1 bins
    # Bandwidths: df = f[i+1] - f[i]
    # Centers: f_mid = (f[i+1] + f[i]) / 2
    dfs = diff(freqs)
    central_freqs = (freqs[1:end-1] + freqs[2:end]) ./ 2.0
    
    # 2. Sample the density from the continuous model at central frequencies
    densities = [ContinuousSpectrums.get_density(shape, f) for f in central_freqs]

    # 4. Compute moments with discretization error
    m₀ = sum(densities .* dfs)

    # 5. Normalization factor (Energy Preservation)
    norm_factor = (shape.Hs / (4.0 * sqrt(m₀)))^2
    
    if mess
        comp_Hs = 4.0 * sqrt(m₀)
        m₁ = sum(central_freqs .* densities .* dfs)
        comp_fs = m₁ / m₀
        println("--- Discrete Spectrum ---")
        println("Spectrum:     ", typeof(shape))
        println("Sampling:     ", typeof(strat))
        println("Frequency range: [", round(fmin, digits=3), ", ", round(fmax, digits=3), "] Hz")
        println("Samples:      ", nf)
        println("Bins (nf):    ", nf-1)
        println("Discrete Hs:  ", round(comp_Hs, digits=3), " m")
        println("Target Hs:    ", round(shape.Hs, digits=3), " m")
        println("Discrete Tp:  ", round(1.0 / comp_fs, digits=3), " s")
        println("Target Tp:    ", round(shape.Tp, digits=3), " s")
        println("Correction:   ", round(norm_factor, digits=3))
        println("----------------------------")
    end

    return DiscreteSpectrum(shape, strat, domain, fmin, fmax, nf, nf-1, norm_factor)
end


# --- GETTERS (Lazy Data Generation) ---

# --- Frequency sampling getters ---

"""
    get_frequencies(spec::DiscreteSpectrum)

Returns the discrete frequency bins and their widths.
"""
function get_frequencies(spec::DiscreteSpectrum)
    return FrequencySampling.generate_grid(spec.sampling, spec.domain, spec.spectrum, spec.fmin, spec.fmax, spec.nf)
end

function get_frequencies(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    all_freqs = get_frequencies(spec)
    return all_freqs[r]
end

function get_frequencies(spec::DiscreteSpectrum, f_range::AbstractRange{<:Real})
    return get_frequencies(spec, get_spectral_index(spec, f_range))
end

function get_frequency(spec::DiscreteSpectrum, f_target::Real)
    # Early return for empty or invalid frequency
    (f_target < spec.fmin || f_target > spec.fmax) && return Float64[]
    return get_frequency(spec, get_spectral_index(spec, f_target))
end

function get_frequency(spec::DiscreteSpectrum, idx::Int)
    all_freqs = get_frequencies(spec)
    return (idx < 1 || idx > spec.nf) ? Float64[] : all_freqs[idx]
end

# -------------------------------------

# --- Central frequency getters ---

"""
    get_central_frequencies(spec::DiscreteSpectrum)

Returns the central frequencies of each discrete bin.
"""
function get_central_frequencies(spec::DiscreteSpectrum)
    freqs = get_frequencies(spec)
    return (freqs[1:end-1] + freqs[2:end]) ./ 2.0
end

function get_central_frequencies(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    central_freqs = get_central_frequencies(spec)
    return central_freqs[r]
end

function get_central_frequency(spec::DiscreteSpectrum, idx::Int)
    central_freqs = get_central_frequencies(spec)
    return (idx < 1 || idx > spec.nbands) ? Float64[] : central_freqs[idx]
end


# -----------------------------------

# --- Spectral index getters ---
"""
    get_spectral_index(spec::DiscreteSpectrum, f_target::Real)

Returns the index (1 to nf-1) of the bin containing the frequency `f_target`.
"""
function get_spectral_index(spec::DiscreteSpectrum, f_target::Real)
    # Early return for empty or invalid frequency
    (f_target < spec.fmin || f_target > spec.fmax) && return Float64[]

    # 2. searchsortedlast finds the highest index i such that edges[i] <= f_target
    idx = searchsortedlast(get_frequencies(spec), f_target)

    # 3. If f_target is exactly the last edge, it technically falls outside 
    # the last bin (since there are only nf-1 bins). 
    # We usually cap it to the last bin index.
    if idx >= spec.nf
        return spec.nbands
    end

    return idx
end


"""
    get_spectral_index(spec::DiscreteSpectrum, f_range::AbstractRange{<:Real})

Returns a UnitRange of bin indexes that cover the provided frequency range.
Useful for extracting a "slice" of the spectrum (e.g., around the peak).
"""
function get_spectral_index(spec::DiscreteSpectrum, f_range::AbstractRange{<:Real})
    f_start = first(f_range)
    f_end   = last(f_range)
    
    # Find start and end bin
    idx_start = get_spectral_index(spec, f_start)
    idx_end   = get_spectral_index(spec, f_end)
    
    # Handle out-of-bounds cases
    if idx_start == 0 && f_start < get_frequencies(spec)[1]
        idx_start = 1
    end
    
    if idx_end == 0 && f_end > get_frequencies(spec)[end]
        idx_end = spec.nf - 1
    end

    # Return as a UnitRange for easy slicing (e.g., spec.f_centers[r])
    return idx_start:idx_end
end

# -----------------------------------

# --- Bandwidth getters ---

"""
    get_bandwidths(spec::DiscreteSpectrum)

Returns the bandwidths of each discrete bin.
"""
function get_bandwidths(spec::DiscreteSpectrum) 
    return diff(get_frequencies(spec))
end

function get_bandwidths(freqs::Vector{Float64}) 
    return diff(freqs)
end

function get_bandwidths(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    dfs = get_bandwidths(spec)
    return dfs[r]
end

function get_bandwidth(spec::DiscreteSpectrum, idx::Int)
    dfs = get_bandwidths(spec)
    return (idx < 1 || idx > spec.nbands) ? Float64[] : dfs[idx]
end

# -------------------------

# --- SPECTRAL DENSITY EVALUATION FUNCTIONS ---

""" 
    get_density(spec::DiscreteSpectrum, idx::Int)

Returns corrected spectral density S(fᵢ) for the specified bin index.
"""
function get_density(spec::DiscreteSpectrum, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]
    # Get the central frequency for the bin
    fᵢ = get_central_frequency(spec, idx)
    # Calculate density on the requested bin
    return ContinuousSpectrums.get_density(spec.spectrum, fᵢ) * spec.norm_factor
end

"""
    get_densities(spec::DiscreteSpectrum)

Returns corrected spectral densities S(fᵢ) for all bins.
"""
function get_densities(spec::DiscreteSpectrum)
    return [ContinuousSpectrums.get_density(spec.spectrum, f) * spec.norm_factor for f in get_central_frequencies(spec)]
end

"""
    get_densities(spec::DiscreteSpectrum, idx0::Int=1, idx1::Int=length(spec.fᵢ))

Returns corrected spectral densities S(fᵢ) for the specified bin range.
"""
function get_densities(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    # Safety: Clamp indexes to valid 1-based range
    start_idx = max(1, first(r))
    end_idx   = min(last(r), spec.nbands)

    # Early return for empty or invalid range
    (start_idx > end_idx) && return Float64[]

    # Calculate densities only for the requested range
    return [ContinuousSpectrums.get_density(spec.spectrum, f) * spec.norm_factor for f in get_central_frequencies(spec, start_idx:end_idx)]
end



# function get_densities(spec::DiscreteSpectrum)
#     # Compute densities at bin edges
#     S_edges = [ContinuousSpectrums.get_density(spec.spectrum, f) * spec.norm_factor for f in get_frequencies(spec)]
#     # Trapezoidal rule
#     return (S_edges[1:end-1] .+ S_edges[2:end]) ./ 2.0
# end

# ---------------------------------------------

# --- SPECTRAL AMPLITUDE EVALUATION FUNCTIONS ---

""" 
    get_amplitude(spec::DiscreteSpectrum, idx::Int)

Returns corrected wave amplitude A(fᵢ) for the specified bin index.
"""
function get_amplitude(spec::DiscreteSpectrum, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]

    # Calculate amplitude only for the requested index
    return sqrt(2.0 * get_density(spec, idx) * get_bandwidth(spec, idx))
end


function get_amplitudes(spec::DiscreteSpectrum)
    return sqrt.(2.0 .* get_densities(spec) .* get_bandwidths(spec))
end

"""
    get_amplitudes(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})

Returns corrected wave amplitudes A(fᵢ) for the specified bin range.
"""
function get_amplitudes(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    # Safety: Clamp indexes
    start_idx = max(1, first(r))
    end_idx   = min(last(r), spec.nbands)

    (start_idx > end_idx) && return Float64[]

    # Calculate amplitudes only for the requested range
    return sqrt.(2.0 .* get_densities(spec, start_idx:end_idx) .* get_bandwidths(spec, start_idx:end_idx))
end

# -----------------------------------------------

"""
    get_moment(spec::DiscreteSpectrum, n::Int)

Calculates the n-th corrected discrete spectral moment: mₙ = Σ (fⁿ * S * df)
"""
function get_moment(spec::DiscreteSpectrum, n::Int)
    vals = [ (f^n) * S for (f, S) in zip(get_central_frequencies(spec), get_densities(spec)) ]
    return sum(vals .* get_bandwidths(spec))
end


"""
    get_Hs(spec::DiscreteSpectrum)

Returns the corrected significant wave height Hs = 4√m₀ computed from the discrete bins.
"""
function get_Hs(spec::DiscreteSpectrum)
    m₀ = get_moment(spec, 0)
    return 4.0 * sqrt(m₀)
end


"""
    discrete_cdf(spec::DiscreteSpectrum, f_target::Real)

Returns the cumulative energy (variance) contained in the spectrum up to 
the frequency `f_target`. 
- If f_target < f_min, returns 0.0.
- If f_target > f_max, returns the total variance (m0).
"""
function discrete_cdf(spec::DiscreteSpectrum, f_target::Real)
    # 1. Boundary Guards
    if f_target < spec.fmin
        return 0.0
    elseif f_target >= spec.fmax
        return get_moment(spec, 0)
    end

    # 3. Binary search to find the index of the largest frequency <= f_target
    # searchsortedlast is a built-in Julia function that returns the index
    idx = get_spectral_index(spec, f_target)
    if idx > 0
        idx = idx - 1 # Adjust to get the last bin fully below f_target
    end
    
    # 4. Return the energy up to that bin
    # Note: For more precision, one could linearly interpolate between 
    # E_cum[idx] and E_cum[idx+1], but in a discrete spectrum context, 
    # returning the value at the bin edge is standard.
    return cumsum(get_densities(spec, 0:idx) .* get_bandwidths(spec, 0:idx))[end]
end

end # module
