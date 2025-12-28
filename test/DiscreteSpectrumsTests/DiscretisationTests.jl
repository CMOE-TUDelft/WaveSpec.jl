using Test
using WaveSpec

@testset "Spectral Spreading - JONSWAP" begin
    # Parameters
    Hs_target = 3.0
    Tp = 10.0
    
    # 1. Create a JONSWAP instance
    # This should trigger our Ag normalization via Gauss Quadrature
    spec = JONSWAP(Hs_target, Tp)
    
    # 2. Test: Check if Ag was calculated
    @test spec.Ag > 0
    
    # 3. Test: Integration of the density
    # The area under S(f) should be (Hs/4)^2
    # m0 = ∫ S(f) df
    f_samples = range(1e-4, 2.0, length=1000)
    df = step(f_samples)
    s_values = [get_density(spec, f) for f in f_samples]
    
    m0 = sum(s_values .* df)
    Hs_calc = 4 * sqrt(m0)
    
    @info "Target Hs: $Hs_target, Calculated Hs: $Hs_calc"
    @test Hs_calc ≈ Hs_target rtol=1e-2
end

@testset "Discrete Spectrum Bins" begin
    spec_base = JONSWAP(3.0, 10.0)
    # Create a discrete container with 100 frequencies
    # Using the bin-edge logic we implemented (Midpoint/Central)
    ds = DiscreteSpectrum(spec_base, f_min=0.02, f_max=0.5, nf=101)
    
    # Check that we have nf-1 bins
    @test length(get_bin_centers(ds)) == 100
    
    # Ensure bandwidths sum to the total range
    @test sum(get_bandwidths(ds)) ≈ (0.5 - 0.02)
end