module BretschneiderTests

using Test
using WaveSpec.ContinuousSpectrums

@testset "Bretschneider Spectrum Specific Tests" begin

    # 1. Setup a standard sea state
    Hs = 2.5
    Tp = 8.0
    fp = 1.0 / Tp
    spec = Bretschneider(Hs, Tp)

    @testset "Mathematical Form & Peak Density" begin
        # For Bretschneider, the peak density S(fp) has a fixed theoretical value:
        # S(fp) = (5/4) * (Hs/4)^2 / (fp * exp(1.25)) * 4
        # Simplified: S(fp) ≈ 0.3125 * (Hs^2 / fp) * exp(-1.25)
        
        val_at_peak = get_density(spec, fp)
        theoretical_peak = 0.3125 * (Hs^2 / fp) * exp(-1.25)
        
        @test val_at_peak ≈ theoretical_peak rtol=1e-5
        
        # Verify it is indeed a peak (local maximum)
        @test val_at_peak > get_density(spec, fp * 0.99)
        @test val_at_peak > get_density(spec, fp * 1.01)
    end

    @testset "Tail Decay (f⁻⁵ Law)" begin
        # Bretschneider must follow S(f) ∝ f⁻⁵ at high frequencies
        f1 = 10.0 * fp
        f2 = 11.0 * fp
        
        s1 = get_density(spec, f1)
        s2 = get_density(spec, f2)
        
        @test s2 / s1 ≈ (f2 / f1)^-5 rtol=1e-4
    end

    @testset "Low Frequency Limit" begin
        # Energy should vanish extremely quickly below peak due to exp(-f⁻⁴)
        @test get_density(spec, 0.1 * fp) ≈ 0.0 atol=1e-10
        @test get_density(spec, 0.0) == 0.0
    end

    @testset "Self-Normalizing Property" begin
        # Unlike JONSWAP, Bretschneider doesn't need an Ag factor integration.
        # We test if the raw formula correctly returns the target Hs.
        # This tests the '0.3125' (which is 1.25/4) constant integrity.
        
        m0_num = integrate(spec) # Using your generic integrate function
        Hs_num = 4 * sqrt(m0_num)
        
        @test Hs_num ≈ Hs rtol=1e-3
    end

end

end

