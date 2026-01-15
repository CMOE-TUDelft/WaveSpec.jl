module DiscretisationTests

using Test
using WaveSpec.ContinuousSpectrums
using WaveSpec.FrequencySampling
using WaveSpec.SpectralSpreading
using WaveSpec.Signal
using Statistics


@testset "Discrete Energy Conservation (Amplitude Test)" begin
    Hs_target = 3.0
    spec = JONSWAP(Hs_target, 10.0)
    
    # Try different sampling strategies
    strategies = [UniformSampling(), LogSampling(), ChebyshevSampling()]
    domains = [Frequency, Energy] 

    for strat in strategies, dom in domains
        # Discretize
        ds = DiscreteSpectrum(spec, strat, get_fmin(spec), get_fmax(spec), 200; domain=dom, mess=false)
        
        # 1. Extract amplitudes
        amplitudes = get_amplitudes(ds)
        
        # 2. Compute total discrete variance (m0)
        # m0 = Σ (1/2 * Ai^2)
        m0_discrete = sum(0.5 .* amplitudes .^ 2)
        
        # 3. Compare to theoretical (Hs/4)^2
        m0_theoretical = (Hs_target / 4.0)^2
        
        # Check if they match within a tight tolerance
        @test m0_discrete ≈ m0_theoretical rtol=1e-3
        
        # 4. Cross-check Hs
        Hs_discrete = 4.0 * sqrt(m0_discrete)
        @test Hs_discrete ≈ Hs_target rtol=1e-3
    end
end


@testset "Uniform Energy Domain Sampling Check" begin
    # 1. Setup
    Hs = 3.0
    spec = JONSWAP(Hs, 10.0)
    nf = 50
    
    # 2. Discretize using the Energy Domain marker
    # We use UniformSampling here to get 'Equal Energy' steps
    ds = DiscreteSpectrum(spec, UniformSampling(), get_fmin(spec), get_fmax(spec), nf; 
                          domain=FrequencySampling.Energy, mess=false)
    
    # 3. Integrate energy per bin
    # All bins should have the same energy stored in them for the uniform energy sampling 
    n_substeps = 200
    f_edges = get_frequencies(ds)
    energies = zeros(ds.nbands)
    for i in 1:ds.nbands
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
            S0 = ContinuousSpectrums.get_density(spec, f_micro[j])
            S1 = ContinuousSpectrums.get_density(spec, f_micro[j+1])
            
            # Area of a small trapezoid: (average height) * width
            bin_energy += 0.5 * (S0 + S1) * df
        end
        # Save energy
        energies[i] = bin_energy
    end

    # We check the relative standard deviation (Coefficient of Variation)
    avg_energ = mean(energies)
    std_energ = std(energies)
    coeff_of_variation = std_energ / avg_energ
    
    @info "Energy Consistency" Mean=avg_energ Std=std_energ CV=coeff_of_variation
    
    # The energies should be extremely consistent
    @test coeff_of_variation < 1e-3
    
    # 5. Check the theoretical value
    # m0 = (Hs/4)^2. Each bin should have m0 / (nf-1) variance.
    expected_energy = (Hs/4.0)^2 / ds.nbands
    
    @test all(isapprox.(energies, expected_energy, rtol=1e-2))
end


@testset "Spectrum Reconstruction (IFT/FFT)" begin
    # --- 1. Discrete Spectrum ---
    Hs_target = 3.0
    Tp_target = 10.0
    spec = JONSWAP(Hs_target, Tp_target)

    # Continuous Spectrum 
    f_cont = range(get_fmin(spec), get_fmax(spec), length=1000)
    S_cont = [ContinuousSpectrums.get_density(spec, f) for f in f_cont]
    
    # Test with Equal Energy sampling to prove the strategy works
    ds = DiscreteSpectrum(spec, UniformSampling(), get_fmin(spec), get_fmax(spec), 750; 
                          domain=FrequencyDomain(), mess=false)
    
    # --- 2. Signal Synthesis ---
    # Sampling parameters
    N = 2^14                        # Even number (power of 2) of samplings in time (choose high number for long time window)
    h = 2.0                         # water depth
    fs = 2.0 *get_fmax(ds.spectrum) # Nyquist frequency (avoid aliasing)
    t, η_signal = generate_signal(ds, h, N, fs=fs)
    
    # --- 3. Statistical Validation ---
    Hs_sim = 4.0 * std(η_signal)
    @test Hs_sim ≈ Hs_target rtol=0.02  # Allow 2% error for stochastic sampling
    
    # --- 4. Spectral Validation (Reconstruction) ---
    f_axis, fAmp, psd = get_single_sided_spectrum(η_signal, fs)
    
    # Test 4.1: Check the peak frequency
    f_peak_sim = f_axis[argmax(psd)]
    f_target = 1.0 / Tp_target
    @test f_peak_sim ≈ f_target atol=0.01
    
    # Test 4.2: The integral of the signal PSD should match the target variance
    @test sum(psd .* fs/N) ≈ Hs_target^2 / 16 rtol=1e-2

    # Test 4.3: Compare Integrated Energy in Bands.
    # For example, check if the energy between $0.8 f_p$ and $1.2 f_p$ is the same in both the model and the simulation:

    # Energy in the peak region
    mask_sim = (f_axis .> 0.08) .& (f_axis .< 0.12)
    energy_sim = sum(psd[mask_sim]) * (f_axis[2] - f_axis[1])

    mask_target = (f_cont .> 0.08) .& (f_cont .< 0.12)
    energy_target = sum(S_cont[mask_target]) * (f_cont[2] - f_cont[1])

    @test energy_sim ≈ energy_target rtol=0.05

    # Test 4.4: Check if the recovered peak density matches the target density
    S_max_target = ContinuousSpectrums.get_density(spec, f_target)
    # Smooth the signal spectrum by averaging over time windows
    f_avg, avg_psd = averaged_psd(f_axis, psd)
    S_max_sim = maximum(avg_psd)
    
    @info "Reconstruction Results" H_target Hs_sim f_peak_target f_peak_sim S_max_target S_max_sim
    
    @test S_max_sim ≈ S_max_target rtol=0.15 
end

end # module