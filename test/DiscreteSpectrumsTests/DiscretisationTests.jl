module DiscretisationTests

using Test
using WaveSpec.ContinuousSpectrums
using WaveSpec.FrequencySampling
using WaveSpec.SpectralSpreading
using Statistics
using FFTW


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
    
    # Test with Equal Energy sampling to prove the strategy works
    ds = DiscreteSpectrum(spec, UniformSampling(), get_fmin(spec), get_fmax(spec), 750; 
                          domain=FrequencyDomain(), mess=false)
    
    # --- 2. Signal Synthesis ---
    # Sampling parameters
    N = 2^14                    # Even number (power of 2) of samplings in time (choose high number for long time window)
    fs = 2.0 * get_fmax(spec)   # Nyquist frequency (the window sampling frequency fs must be MORE than double the signal maximum frequency -> AVOID ALIASING)
    dt = 1.0 / fs
    t = collect(0:dt:(N-1)*dt)
    
    amplitudes = SpectralSpreading.get_amplitudes(ds)
    freqs = SpectralSpreading.get_central_frequencies(ds)
    phases = 2π .* rand(ds.nf)
    
    # η(t) = Σ A_i * cos(2π f_i t + ϕ_i)
    eta = zeros(length(t))
    for (A, f, ϕ) in zip(amplitudes, freqs, phases)
        eta .+= A .* cos.(2π .* f .* t .+ ϕ)
    end
    
    # --- 3. Statistical Validation ---
    Hs_sim = 4.0 * std(eta)
    @test Hs_sim ≈ Hs_target rtol=0.02  # Allow 2% error for stochastic sampling
    
    # --- 4. Spectral Validation (Reconstruction) ---
    # Perform the signal's Fast Fourier Transform
    fft_res = fft(eta)
    # Power Spectral Density (PSD)
    # Factor 2.0 because we use only positive frequencies
    psd_sim = (2.0 * dt / N) .* abs.(fft_res[1:Int(N/2)]) .^ 2
    f_axis = (0:Int(N/2)-1) ./ last(t)
    
    # Check the peak frequency
    f_peak_sim = f_axis[argmax(psd_sim)]
    @test f_peak_sim ≈ (1/Tp_target) atol=0.01
    
    # Check if the recovered peak density matches the target density
    # We compare the peak of the simulated PSD to the peak of the JONSWAP model
    S_max_target = ContinuousSpectrums.get_density(spec, 1/Tp_target)
    S_max_sim = maximum(psd_sim)
    
    @info "Reconstruction Results" Hs_sim f_peak_sim S_max_sim S_max_target
    
    # Note: PSD from FFT can be noisy. In a real test, 
    # we might use Welch's method (averaging) for a smoother comparison.
    @test S_max_sim ≈ S_max_target rtol=0.15 
end

end # module