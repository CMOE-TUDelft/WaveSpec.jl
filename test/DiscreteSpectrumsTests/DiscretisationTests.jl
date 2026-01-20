module DiscretisationTests

using Test
using WaveSpec.ContinuousSpectrums
using WaveSpec.SpectralSampling
using WaveSpec.SpectralSpreading
using WaveSpec.Signal
using Statistics


"""
    prepare_test_models()
Prepares a list of implemented continuous spectral models, with tuples: (model_name::String, model_instance::AbstractSpectrum)
"""
function prepare_test_models()
    # Target Parameters
    Hs_target = 3.0
    Tp_target = 2.0
    h_target  = 30.0
    U10_target = 10.0

    Hs_target_2 = 1.0
    Tp_target_2 = 4.0
    λ1 = 2
    λ2 = 3

    # Prepare Models in list of tuples (name, structure)
    models = [("JONSWAP", JONSWAP(Hs_target, Tp_target)), 
              ("TMA", TMA(Hs_target, Tp_target, h_target)),
              ("Donelan", Donelan(Hs_target, Tp_target, U10_target)), 
              ("Bretschneider", Bretschneider(Hs_target, Tp_target)),
              ("OchiHubble", OchiHubble(Hs_target, Tp_target, λ1, Hs_target_2, Tp_target_2, λ2))]

    return models
end


function test_energy_conservation_amplitude(spec::AbstractSpectrum)
    # Extract target parameters
    Hs_target = ContinuousSpectrums.get_Hs(spec)

    # Try different sampling strategies
    strategies = [UniformSampling(), LogSampling(), ChebyshevSampling()]
    domains = [Frequency, Energy] 

    for strat in strategies, dom in domains
        # Discretize
        ds = DiscreteSpectralSpreading(spec, strat, get_fmin(spec), get_fmax(spec), 200; domain=dom, mess=false)
        
        # 1. Extract amplitudes
        amplitudes = get_amplitudes(ds)
        
        # 2. Compute total discrete variance (m0)
        # m0 = Σ (1/2 * Ai^2)
        m0_discrete = sum(0.5 .* amplitudes .^ 2)
        
        # 3. Compare to theoretical (Hs/4)^2
        m0_theoretical = (Hs_target / 4.0)^2
        
        # Check if they match within a tight tolerance
        #@info "1srt order moment" Theoretical=m0_theoretical Discrete=m0_discrete 
        @test m0_discrete ≈ m0_theoretical rtol=1e-3
        
        # 4. Cross-check Hs
        Hs_discrete = 4.0 * sqrt(m0_discrete)
        #@info "Significant Wave Height" Target=Hs_target Discrete=Hs_discrete
        @test Hs_discrete ≈ Hs_target rtol=1e-3
    end
end



function test_uniform_energy_sampling(spec::AbstractSpectrum)
    Hs_target = ContinuousSpectrums.get_Hs(spec)

    # Discretize using the Energy Domain marker
    nf = 20
    ds = DiscreteSpectralSpreading(spec, UniformSampling(), get_fmin(spec), get_fmax(spec), nf; 
                          domain=SpectralSampling.Energy, mess=false)
    
    # Integrate energy per bin
    f_edges = get_frequencies(ds)
    energies = zeros(ds.nbands)
    for i in 1:ds.nbands
        # 1. Define the bin boundaries
        f_start = f_edges[i]
        f_end   = f_edges[i+1]
        
        # 2. Integrate energy in the bin
        energies[i] = integrate(spec, f_start, f_end)
    end

    # We check the relative standard deviation (Coefficient of Variation)
    avg_energ = mean(energies)
    std_energ = std(energies)
    coeff_of_variation = std_energ / avg_energ
    
    #@info "Energy Consistency" Mean=avg_energ Std=std_energ CV=coeff_of_variation
    
    # The energies should be extremely consistent
    @test coeff_of_variation < 1e-3
    
    # 5. Check the theoretical value
    # m0 = (Hs/4)^2. Each bin should have m0 / (nf-1) variance.
    expected_energy = (Hs_target/4.0)^2 / ds.nbands
    
    @test all(isapprox.(energies, expected_energy, rtol=1e-2))
end



function test_signal_reconstruction(spec::AbstractSpectrum)
    # Extract target parameters
    Hs_target = ContinuousSpectrums.get_Hs(spec)
    fp_target = 1.0 / ContinuousSpectrums.get_Tp(spec)

    # Continuous Spectrum 
    f_cont = range(get_fmin(spec), get_fmax(spec), length=1000)
    S_cont = [ContinuousSpectrums.get_density(spec, f) for f in f_cont]
    
    # Test with Equal Energy sampling to prove the strategy works
    nf = 2^13
    ds = DiscreteSpectralSpreading(spec, UniformSampling(), get_fmin(spec), get_fmax(spec), nf; 
                          domain=FrequencyDomain(), mess=false)
    
    # --- 2. Signal Synthesis ---
    # Sampling parameters
    N = 2^14   # N = 2 * nf         # Even number (power of 2) of samplings in time (choose high number for long time window)
    h = 2.0                         # water depth
    fs = 2.0 *get_fmax(ds.spectrum) # Nyquist frequency (avoid aliasing)
    t, η_signal = generate_signal(ds, h, N, fs=fs)
    
    # --- 3. Statistical Validation ---
    Hs_sim = 4.0 * std(η_signal)
    @test Hs_sim ≈ Hs_target rtol=0.02  # Allow 2% error for stochastic sampling
    
    # --- 4. Spectral Validation (Reconstruction) ---
    f_axis, fAmp, psd = get_single_sided_spectrum(η_signal, fs)
    # Smooth the signal spectrum by averaging over time windows
    f_avg, avg_psd = averaged_psd(f_axis, psd)
    
    # Test 4.1: Check the peak frequency
    fp_sim = f_avg[argmax(avg_psd)]
    @test fp_sim ≈ fp_target atol=0.01
    
    # Test 4.2: The integral of the signal PSD should match the target variance
    @test sum(psd .* fs/N) ≈ Hs_target^2 / 16 rtol=1e-2

    # Test 4.3: Compare Integrated Energy in Bands.
    mask_sim = (f_axis .> 0.8*fp_target) .& (f_axis .< 1.2*fp_target)
    energy_sim = sum(psd[mask_sim]) * (f_axis[2] - f_axis[1])

    mask_target = (f_cont .> 0.8*fp_target) .& (f_cont .< 1.2*fp_target)
    energy_target = sum(S_cont[mask_target]) * (f_cont[2] - f_cont[1])

    @test energy_sim ≈ energy_target rtol=0.05

    # Test 4.4: Check if the recovered peak density matches the target density
    S_max_target = ContinuousSpectrums.get_density(spec, fp_target)
    S_max_sim = maximum(avg_psd)
    
    #@info "Reconstruction Results" Hs_target Hs_sim fp_target fp_sim S_max_target S_max_sim
    
    @test S_max_sim ≈ S_max_target rtol=0.15 
end



@testset "Discrete Energy Conservation (Amplitude Test)" begin
    # Prepare models
    models = prepare_test_models()
    # Test models
    for (name, model) in models
        @testset "Spectrum: $name" begin    
            test_energy_conservation_amplitude(model)
        end
    end
end

@testset "Uniform Energy Domain Sampling Check" begin
    # Prepare models
    models = prepare_test_models()
    # Test models
    for (name, model) in models
        @testset "Spectrum: $name" begin    
            test_uniform_energy_sampling(model)
        end
    end
end

@testset "Spectrum Reconstruction (IFT/FFT)" begin
    # Prepare models
    models = prepare_test_models()
    # Test models
    for (name, model) in models
        @testset "Spectrum: $name" begin    
            test_signal_reconstruction(model)
        end
    end
end

end # module