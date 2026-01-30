module TMATests

using Test
using WaveSpec.ContinuousSpectrums

@testset "TMA Spectrum Specific Tests" begin 

    # 1. Setup Parameters
    Hs = 2.0
    Tp = 8.0
    h_shallow = 5.0   # Shallow water (depth < wavelength)
    h_deep = 500.0    # Deep water (depth >> wavelength)

    spec_shallow = TMA(Hs, Tp, h_shallow)
    spec_deep = TMA(Hs, Tp, h_deep)

    @testset "Kitaigorodskii Phi Factor Limits" begin
        # Test the piecewise function get_phi(f, h)
        # 1. Deep water / High frequency: Phi should approach 1.0
        @test get_phi(spec_deep, 10.0) ≈ 1.0
        
        # 2. Very shallow / Low frequency: Phi should be < 1.0
        # For f=0.1, h=5, ωh ≈ 0.628 * sqrt(5/9.81) ≈ 0.44 -> returns 0.5 * ωh^2
        @test get_phi(spec_shallow, 0.1) < 1.0
        @test get_phi(spec_shallow, 0.0) == 0.0
    end

    @testset "Depth Sensitivity (Shallow vs Deep)" begin
        # For the same Hs and Tp, the TMA spectrum in shallow water 
        # needs a HIGHER Ag to reach the same total energy because 
        # the Phi factor is removing energy from the tail.
        @test spec_shallow.Ag_tma > spec_deep.Ag_tma
        
        # In very deep water, TMA should converge toward JONSWAP
        # The Ag_tma should be very close to the underlying JONSWAP Ag
        @test spec_deep.Ag_tma ≈ spec_deep.js.Ag rtol=1e-3

        # For shallow water, Ag_tma should be significantly larger
        # than the JONSWAP Ag to compensate for tail loss
        @test spec_shallow.Ag_tma > 1.1 * spec_shallow.js.Ag
    end

    @testset "High Frequency Tail Attenuation" begin
        # In shallow water, the f^-5 tail of JONSWAP is modified.
        # Check that high frequency density in shallow water is 
        # lower than deep water before Ag normalization is applied.
        f_high = 2.0 * (1.0 / Tp)
        
        # We check the raw impact of Phi
        phi_shallow = get_phi(spec_shallow, f_high)
        phi_deep = get_phi(spec_deep, f_high)
        
        @test phi_shallow < phi_deep
        @test phi_deep ≈ 1.0
    end

end

end