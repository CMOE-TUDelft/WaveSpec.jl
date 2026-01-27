module DonelanTests

using Test
using WaveSpec.PhysicalConstants: g
using WaveSpec.ContinuousSpectrums

@testset "Donelan Spectrum Specific Tests" begin

    # 1. Setup Parameters
    Hs = 2.0
    Tp = 6.0
    U10 = 15.0 # Strong wind
    spec = Donelan(Hs, Tp, U10)

    @testset "Wind Speed & Wave Age Logic" begin
        # Test γ dependency on U10
        # Phase speed cp = g / (2π * fp)
        fp = 1.0 / Tp
        cp = g / (2π * fp)
        U_cp = U10 / cp
        
        # If U_cp > 1, γ should be > 1.7
        if U_cp > 1.0
            @test spec.γ > 1.7
        else
            @test spec.γ == 1.7
        end

        # Test limits: check if γ stays within the 1.7 - 7.0 bounds
        spec_slow = Donelan(Hs, Tp, 5.0)  # Very low wind
        spec_fast = Donelan(Hs, Tp, 50.0) # Hurricane force
        @test spec_slow.γ >= 1.7
        @test spec_fast.γ <= 7.0
    end

    @testset "Tail Decay (f⁻⁴ Law)" begin
        # This is the most critical test for Donelan.
        # At high frequencies, S(f) ∝ f⁻⁴
        f1 = 10.0 * spec.fp
        f2 = 11.0 * spec.fp
        
        s1 = get_density(spec, f1)
        s2 = get_density(spec, f2)
        
        @test s2 / s1 ≈ (f2 / f1)^-4 rtol=1e-2
    end

    @testset "Shape Parameters (σ₁ and σ₂)" begin
        # Donelan uses 0.08 and 0.10 (slightly wider than JONSWAP)
        @test spec.σ₁ == 0.08
        @test spec.σ₂ == 0.10
        
        # Verify asymmetry: right side (σ₂) is broader than left side (σ₁)
        fp = spec.fp
        @test get_density(spec, fp + 0.1) > get_density(spec, fp - 0.1)
    end

end

end