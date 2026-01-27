module JONSWAPTests

using Test
using WaveSpec.ContinuousSpectrums

@testset "JONSWAP Spectrum Specific Tests" begin

    Hs, Tp = 3.0, 10.0

    
    # A unique feature of JONSWAP is that the density at the peak f = f_p is exactly γ times 
    # higher than the "base" PM-style spectrum at that same frequency (provided the normalization 
    # Ag is handled). You can verify that as γ increases, the peak density increases proportionally.
    @testset "JONSWAP-Specific: γ Peak Enhancement" begin
        # Compare a PM-limit (γ=1) to a standard JONSWAP (γ=3.3)
        spec_gamma1 = JONSWAP(Hs, Tp, 1.0)
        spec_gamma3 = JONSWAP(Hs, Tp, 3.3)
        
        # At the peak fr = 1.0, the enhancement part is γ^exp(0) = γ^1
        # Note: Ag also changes, so we check the ratio of the peak densities
        peak_ratio = get_density(spec_gamma3, 1/Tp) / get_density(spec_gamma1, 1/Tp)
        
        @test peak_ratio > 1.0
        @test spec_gamma3.γ == 3.3
    end

    
    # JONSWAP is defined by two different widths: σ = 0.07 (left of peak) and σ = 0.09 (right of peak). 
    # This means the spectrum is steeper on the left than the right.
    @testset "JONSWAP-Specific: Spectral Asymmetry" begin
        spec = JONSWAP(Hs, Tp)
        fp = spec.fp
        df = 0.02 # small offset
        
        val_left  = get_density(spec, fp - df)
        val_right = get_density(spec, fp + df)
        
        # Because σ₁ (0.07) < σ₂ (0.09), the decay to the left is sharper.
        # Therefore, the density at the same distance to the right should be higher.
        @test val_right > val_left
    end

    
    # The normalization factor Ag must decrease as γ increases to keep the total energy (m0)
    # constant for the same Hs and Tp.
    @testset "JONSWAP-Specific: Ag Normalization" begin
        # Ag must decrease as γ increases to keep m0 = (Hs/4)² constant
        spec_1 = JONSWAP(Hs, Tp, 1.0)
        spec_5 = JONSWAP(Hs, Tp, 5.0)
        
        @test spec_5.Ag < spec_1.Ag
    end

    @testset "JONSWAP-Specific: Tail Decay" begin
        spec = JONSWAP(Hs, Tp)
        f1 = 5.0 * spec.fp
        f2 = 10.0 * spec.fp
        
        s1 = get_density(spec, f1)
        s2 = get_density(spec, f2)
        
        # S2 / S1 should be roughly (f2 / f1)^-5 = (2)^-5 = 1/32
        observed_ratio = s2 / s1
        expected_ratio = (f2 / f1)^-5
        
        @test observed_ratio ≈ expected_ratio rtol=1e-2
    end

end

end # module