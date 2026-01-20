module SpectralSpreading

using ..SpectralSampling
using ..ContinuousSpectrums
using ..Integration
using RecipesBase

export DiscreteSpectralSpreading
export get_frequency, get_frequencies,
       get_density, get_densities,
       get_amplitude, get_amplitudes 
export get_moments, get_Hs, get_energies

"""
    DiscreteSpectralSpreading
The discrete realization of a frequency spectrum. 
Represents the energy distribution across a finite set of frequency bins.
"""
struct DiscreteSpectralSpreading
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
    DiscreteSpectralSpreading(shape, strategy, fmin, fmax, nf)

Continuous abstract spectrum discretised according to the selected abstract sampling strategy, between fmin and fmax with nf bins.
"""
function DiscreteSpectralSpreading(shape::AbstractSpectrum, strat::AbstractSampling, fmin::Real, fmax::Real, nf::Int; domain::SamplingDomain = Frequency, mess::Bool=true)
    # 1. Generate the nf edges
    freqs = SpectralSampling.generate_grid(strat, domain, shape, fmin, fmax, nf)
    
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
    norm_factor = (ContinuousSpectrums.get_Hs(shape) / (4.0 * sqrt(m₀)))^2
    
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
        println("Target Hs:    ", round(ContinuousSpectrums.get_Hs(shape), digits=3), " m")
        println("Discrete Tp:  ", round(1.0 / comp_fs, digits=3), " s")
        println("Target Tp:    ", round(ContinuousSpectrums.get_Tp(shape), digits=3), " s")
        println("Correction:   ", round(norm_factor, digits=3))
        println("----------------------------")
    end

    return DiscreteSpectralSpreading(shape, strat, domain, fmin, fmax, nf, nf-1, norm_factor)
end


# --- GETTERS (Lazy Data Generation) ---

# --- Frequency sampling getters ---

"""
    get_frequencies(spec::DiscreteSpectralSpreading)

Returns the discrete frequency bins and their widths.
"""
function get_frequencies(spec::DiscreteSpectralSpreading)
    return SpectralSampling.generate_grid(spec.sampling, spec.domain, spec.spectrum, spec.fmin, spec.fmax, spec.nf)
end

function get_frequencies(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
    all_freqs = get_frequencies(spec)
    return all_freqs[r]
end

function get_frequencies(spec::DiscreteSpectralSpreading, f_range::AbstractRange{<:Real})
    return get_frequencies(spec, get_spectral_index(spec, f_range))
end

function get_frequency(spec::DiscreteSpectralSpreading, f_target::Real)
    # Early return for empty or invalid frequency
    (f_target < spec.fmin || f_target > spec.fmax) && return Float64[]
    return get_frequency(spec, get_spectral_index(spec, f_target))
end

function get_frequency(spec::DiscreteSpectralSpreading, idx::Int)
    all_freqs = get_frequencies(spec)
    return (idx < 1 || idx > spec.nf) ? Float64[] : all_freqs[idx]
end

# -------------------------------------

# --- Central frequency getters ---

"""
    get_central_frequencies(spec::DiscreteSpectralSpreading)

Returns the central frequencies of each discrete bin.
"""
function get_central_frequencies(spec::DiscreteSpectralSpreading)
    freqs = get_frequencies(spec)
    return (freqs[1:end-1] + freqs[2:end]) ./ 2.0
end

function get_central_frequencies(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
    central_freqs = get_central_frequencies(spec)
    return central_freqs[r]
end

function get_central_frequency(spec::DiscreteSpectralSpreading, idx::Int)
    central_freqs = get_central_frequencies(spec)
    return (idx < 1 || idx > spec.nbands) ? Float64[] : central_freqs[idx]
end


# -----------------------------------

# --- Spectral index getters ---
"""
    get_spectral_index(spec::DiscreteSpectralSpreading, f_target::Real)

Returns the index (1 to nf-1) of the bin containing the frequency `f_target`.
"""
function get_spectral_index(spec::DiscreteSpectralSpreading, f_target::Real)
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
    get_spectral_index(spec::DiscreteSpectralSpreading, f_range::AbstractRange{<:Real})

Returns a UnitRange of bin indexes that cover the provided frequency range.
Useful for extracting a "slice" of the spectrum (e.g., around the peak).
"""
function get_spectral_index(spec::DiscreteSpectralSpreading, f_range::AbstractRange{<:Real})
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
    get_bandwidths(spec::DiscreteSpectralSpreading)

Returns the bandwidths of each discrete bin.
"""
function get_bandwidths(spec::DiscreteSpectralSpreading) 
    return diff(get_frequencies(spec))
end

function get_bandwidths(freqs::Vector{Float64}) 
    return diff(freqs)
end

function get_bandwidths(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
    dfs = get_bandwidths(spec)
    return dfs[r]
end

function get_bandwidth(spec::DiscreteSpectralSpreading, idx::Int)
    dfs = get_bandwidths(spec)
    return (idx < 1 || idx > spec.nbands) ? Float64[] : dfs[idx]
end

# -------------------------

# --- SPECTRAL DENSITY EVALUATION FUNCTIONS ---

""" 
    get_density(spec::DiscreteSpectralSpreading, idx::Int)

Returns corrected spectral density S(fᵢ) for the specified bin index.
"""
function get_density(spec::DiscreteSpectralSpreading, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]
    # Get the central frequency for the bin
    fᵢ = get_central_frequency(spec, idx)
    # Calculate density on the requested bin
    return ContinuousSpectrums.get_density(spec.spectrum, fᵢ) * spec.norm_factor
end

"""
    get_densities(spec::DiscreteSpectralSpreading)

Returns corrected spectral densities S(fᵢ) for all bins.
"""
function get_densities(spec::DiscreteSpectralSpreading)
    return [ContinuousSpectrums.get_density(spec.spectrum, f) * spec.norm_factor for f in get_central_frequencies(spec)]
end

"""
    get_densities(spec::DiscreteSpectralSpreading, idx0::Int=1, idx1::Int=length(spec.fᵢ))

Returns corrected spectral densities S(fᵢ) for the specified bin range.
"""
function get_densities(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
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
    get_energy(spec::DiscreteSpectralSpreading, idx::Int)

Returns corrected wave energy E(fᵢ) for the specified bin index.
"""
function get_energy(spec::DiscreteSpectralSpreading, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]

    # Calculate amplitude only for the requested index
    return get_density(spec, idx) * get_bandwidth(spec, idx)
end


function get_energies(spec::DiscreteSpectralSpreading)
    return get_densities(spec) .* get_bandwidths(spec)
end

"""
    get_energies(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})

Returns corrected wave energies E(fᵢ) for the specified bins range.
"""
function get_energies(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
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
    get_amplitude(spec::DiscreteSpectralSpreading, idx::Int)

Returns corrected wave amplitude A(fᵢ) for the specified bin index.
"""
function get_amplitude(spec::DiscreteSpectralSpreading, idx::Int)
    # Early return for empty or invalid index
    (idx < 1 || idx > spec.nbands) && return Float64[]

    # Calculate amplitude only for the requested index
    return sqrt(2.0 * get_energy(spec, idx))
end


function get_amplitudes(spec::DiscreteSpectralSpreading)
    return sqrt.(2.0 .* get_energies(spec))
end

"""
    get_amplitudes(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})

Returns corrected wave amplitudes A(fᵢ) for the specified bin range.
"""
function get_amplitudes(spec::DiscreteSpectralSpreading, r::AbstractUnitRange{Int})
    # Safety: Clamp indexes
    start_idx = max(1, first(r))
    end_idx   = min(last(r), spec.nbands)

    (start_idx > end_idx) && return Float64[]

    # Calculate amplitudes only for the requested range
    return sqrt.(2.0 .* get_energies(spec, start_idx:end_idx))
end

# -----------------------------------------------

"""
    get_moment(spec::DiscreteSpectralSpreading, n::Int)

Calculates the n-th corrected discrete spectral moment: mₙ = Σ (fⁿ * S * df)
"""
function get_moment(spec::DiscreteSpectralSpreading, n::Int)
    vals = [ (f^n) * S for (f, S) in zip(get_central_frequencies(spec), get_densities(spec)) ]
    return sum(vals .* get_bandwidths(spec))
end


"""
    get_Hs(spec::DiscreteSpectralSpreading)

Returns the corrected significant wave height Hs = 4√m₀ computed from the discrete bins.
"""
function get_Hs(spec::DiscreteSpectralSpreading)
    m₀ = get_moment(spec, 0)
    return 4.0 * sqrt(m₀)
end

""" 
    get_integrated_energies(spec::DiscreteSpectralSpreading; nfine = 200)

Returns the energies per bin by integrating the continuous spectrum inside
each bin (-> continuous spectrum-wise energy inside discrete-wise bins).
"""
function get_integrated_energies(spec::DiscreteSpectralSpreading; nfine = 200)
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
#   2. Call this function as plot(DiscreteSpectralSpreading)

@recipe function f(spec::DiscreteSpectralSpreading)
    # Set up plot layout
    layout := (1, 2)
    size   := (1000, 400)

    # 1. Continuous Model
    f_range = range(ContinuousSpectrums.get_fmin(spec.spectrum), ContinuousSpectrums.get_fmax(spec.spectrum), length=300)
    S_cont = [ContinuousSpectrums.get_density(spec.spectrum, f) for f in f_range]

    # 2. Get Discrete Data
    f_edges = get_frequencies(spec)                                                      # Edge frequencies
    S_edges = [ContinuousSpectrums.get_density(spec.spectrum, f) for f in f_edges]       # Density at edge frequencies
    f_centers = get_central_frequencies(spec)                                            # Central frequencies
    Δf = get_bandwidths(spec)                                                            # Widths of each bin
    S_centers = [ContinuousSpectrums.get_density(spec.spectrum, f) for f in f_centers]   # Density at central frequencies (S_j)

    collect_f = []
    collect_S = []
    for (f, w) in zip(f_edges, S_edges)
        push!(collect_f, f, f, NaN)
        push!(collect_S, 0, w, NaN)
    end

    # 3. Energy per Bin 
    E_bins = get_integrated_energies(spec)
    
    # --- LEFT SUBPLOT: Spectral Density ---

    # 1a. Continuous spectrum
    @series begin
        subplot := 1
        label := "Continuous Spectrum"
        linecolor := :black
        lw := 2
        f_range, S_cont
    end

    # 1b. Spectral Bins (Bars)
    @series begin
        subplot    := 1
        seriestype := :bar
        label      := "Spectral bins"
        color      := :orange
        alpha      := 0.3
        bar_width  := Δf
        f_centers, S_centers
    end

    # 1c. Vertical lines for bin edges
    @series begin
        subplot    := 1
        seriestype := :vline
        label      := false
        color      := :gray
        alpha      := 0.3
        f_edges
    end

    # 1d. Discrete Samples (Points + Stems)
    @series begin
        subplot   := 1
        label     := false
        linecolor := :red
        lw        := 1
        collect_f, collect_S
    end

    @series begin
        subplot    := 1
        seriestype := :scatter
        label      := "Discrete Samples"
        marker     := :circle
        markersize := 3
        markercolor := :red
        f_edges, S_edges
    end

    # Axis Formatting for Subplot 1
    @series begin
        subplot := 1
        label   := false
        xguide  := "Frequency f [Hz]"
        yguide  := "Spectral Density S(f) [m²s]"
        title   := "Spectral Discretization"
        # This is a dummy series to apply labels to the subplot
        [], []
    end

    # --- SUBPLOT 2: Energy Distribution ---

    @series begin
        subplot    := 2
        seriestype := :bar
        label      := false
        color      := :green
        alpha      := 0.6
        xguide     := "Bin Index"
        yguide     := "Energy E [m²]"
        title      := "Energy Distribution"
        ylims      := (0, maximum(E_bins) * 1.2)
        1:length(E_bins), E_bins
    end

end

end # module
