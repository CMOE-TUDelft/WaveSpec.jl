module OchiHubbleTests

using Test
using WaveSpec.ContinuousSpectrums

@testset "Ochi-Hubble Spectrum Specific Tests" begin

    # 1. Setup a typical bimodal sea state
    # Component 1: Swell (Low frequency, narrow peak λ=1)
    # Component 2: Wind Sea (High frequency, broad peak λ=3)
    Hs1, Tp1, λ1 = 1.5, 14.0, 1.0
    Hs2, Tp2, λ2 = 2.5,  7.0, 3.0

    spec = OchiHubble(Hs1, Tp1, λ1, Hs2, Tp2, λ2)

    @testset "Bimodal Peak Locations" begin
        # Verify that the density is significant at both peak frequencies
        f_swell = 1.0 / Tp1
        f_wind  = 1.0 / Tp2
        
        @test get_density(spec, f_swell) > 0.1
        @test get_density(spec, f_wind) > 0.1
        
        # Verify that a frequency between the peaks has lower density 
        # (assuming peaks are sufficiently separated)
        f_mid = (f_swell + f_wind) / 2.0
        @test get_density(spec, f_mid) < get_density(spec, f_wind)
    end

    @testset "Total Significant Wave Height" begin
        # Ochi-Hubble should store the combined Hs: sqrt(Hs1² + Hs2²)
        expected_total = sqrt(Hs1^2 + Hs2^2)
        @test get_Hs(spec) ≈ expected_total
    end

    @testset "Shape Factor (λ) Influence" begin
        # Higher λ means a narrower, higher peak for the same Hs and Tp
        # Compare two single-component setups
        spec_broad = OchiHubble(2.0, 10.0, 1.0, 0.0, 20.0, 1.0) # λ=1
        spec_sharp = OchiHubble(2.0, 10.0, 5.0, 0.0, 20.0, 1.0) # λ=5
        
        fp = 1.0 / 10.0
        @test get_density(spec_sharp, fp) > get_density(spec_broad, fp)
    end

    @testset "Component Independence" begin
        # If one Hs is set to 0, the total density should equal the other component
        spec_only_swell = OchiHubble(Hs1, Tp1, λ1, 0.0, Tp2, λ2)
        
        f = 1.0 / Tp1
        S_total = get_density(spec_only_swell, f)
        S_swell_only = compute_ochi_component(f, Hs1, 1/Tp1, λ1)
        
        @test S_total ≈ S_swell_only
    end

end

end