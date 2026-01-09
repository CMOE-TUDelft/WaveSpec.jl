module SpectralSpreading

using ..FrequencySampling
using ..ContinuousSpectrums
using ..Integration
using RecipesBase

export DiscreteSpectrum
export get_frequency, get_frequencies,
       get_density, get_densities,
       get_amplitude, get_amplitudes 
export get_moments, get_Hs, get_energies

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

# ---------------------------------------------

# --- SPECTRAL ENERGY EVALUATION FUNCTIONS ---

""" 
    get_energy(spec::DiscreteSpectrum, idx::Int)

Returns corrected wave energy E(fᵢ) for the specified bin index.
"""
function get_energy(spec::DiscreteSpectrum, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]

    # Calculate amplitude only for the requested index
    return get_density(spec, idx) * get_bandwidth(spec, idx)
end


function get_energies(spec::DiscreteSpectrum)
    return get_densities(spec) .* get_bandwidths(spec)
end

"""
    get_energies(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})

Returns corrected wave energies E(fᵢ) for the specified bins range.
"""
function get_energies(spec::DiscreteSpectrum, r::AbstractUnitRange{Int})
    # Safety: Clamp indexes
    start_idx = max(1, first(r))
    end_idx   = min(last(r), spec.nbands)

    (start_idx > end_idx) && return Float64[]

    # Calculate amplitudes only for the requested range
    return get_densities(spec, start_idx:end_idx) .* get_bandwidths(spec, start_idx:end_idx)
end

# --------------------------------------------

# --- SPECTRAL AMPLITUDE EVALUATION FUNCTIONS ---

""" 
    get_amplitude(spec::DiscreteSpectrum, idx::Int)

Returns corrected wave amplitude A(fᵢ) for the specified bin index.
"""
function get_amplitude(spec::DiscreteSpectrum, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]

    # Calculate amplitude only for the requested index
    return sqrt(2.0 * get_energy(spec, idx))
end


function get_amplitudes(spec::DiscreteSpectrum)
    return sqrt.(2.0 .* get_energies(spec))
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
    return sqrt.(2.0 .* get_energies(spec, start_idx:end_idx))
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
    get_integrated_energies(spec::DiscreteSpectrum; nfine = 200)

Returns the energies per bin by integrating the continuous spectrum inside
each bin (-> continuous spectrum-wise energy inside discrete-wise bins).
"""
function get_integrated_energies(spec::DiscreteSpectrum; nfine = 200)
    # Allocate energies array
    e_bins = zeros(spec.nbands)
    # get edge frequencies
    f_edges = get_frequencies(spec)
    # integrate continuous spectrum inside each bin
    for i in 1:spec.nbands
        e_bins[i] = ContinuousSpectrums.integrate(spec.spectrum, f_edges[i], f_edges[i+1]) 
    end
    return e_bins
end



# --- Plots and Visualisation ---
# Plot using light package RecipesBase:
#   1. import Plots in script
#   2. Call this function as plot(DiscreteSpectrum)

@recipe function f(ds::DiscreteSpectrum)
    layout := (1, 2)

    # 1. Continuous Model
    f_range = range(get_fmin(Cspec), get_fmax(Cspec), length=300)
    S_cont = [ContinuousSpectrums.get_density(Cspec, f) for f in f_range]

    # 2. Get Discrete Data
    f_edges = SpectralSpreading.get_frequencies(Dspec)                           # Edge frequencies
    S_edges = [ContinuousSpectrums.get_density(Cspec, f) for f in f_edges]       # Density at edge frequencies
    f_centers = SpectralSpreading.get_central_frequencies(Dspec)                 # Central frequencies
    Δf = diff(f_edges)                                                           # Widths of each bin
    S_centers = [ContinuousSpectrums.get_density(Cspec, f) for f in f_centers]  # Density at central frequencies (S_j)

    collect_f = []
    collect_S = []
    for (f, w) in zip(f_edges, S_edges)
        push!(collect_f, f, f, NaN)
        push!(collect_S, 0, w, NaN)
    end

    # 3. Energy per Bin Calculation
    # We calculate E_j = ∫ S(f) df over the bin limits [f_low, f_high]  using the high-resolution CDF
    E_bins = zeros(Dspec.nbands)
    n_substeps = 100
    for i in 1:Dspec.nbands
        # 1. Define the bin boundaries
        f_start = f_edges[i]
        f_end   = f_edges[i+1]
        
        # 2. Create a micro-grid inside this specific bin
        # We sample n_substeps points to capture the curve's shape
        f_micro = range(f_start, f_end, length=n_substeps)
        df = (f_end - f_start) / (n_substeps - 1)
        
        # 3. Integrate S(f) using the Trapezoidal Rule
        bin_energy = 0.0
        for j in 1:(n_substeps - 1)
            S0 = ContinuousSpectrums.get_density(Cspec, f_micro[j])
            S1 = ContinuousSpectrums.get_density(Cspec, f_micro[j+1])
            
            # Area of a small trapezoid: (average height) * width
            bin_energy += 0.5 * (S0 + S1) * df
        end
        # Save energy
        E_bins[i] = bin_energy
    end
    
    # Left Plot
    @series begin
        subplot := 1
        # ... (your existing density bar code) ...
    end
    
    # Right Plot
    @series begin
        subplot := 2
        seriestype := :bar
        # ... (E_bins logic) ...
    end
end

end # module
