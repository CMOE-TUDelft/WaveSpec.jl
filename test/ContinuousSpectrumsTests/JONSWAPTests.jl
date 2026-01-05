module JONSWAPTests

using Test
using WaveSpec.ContinuousSpectrums
using WaveSpec.Integration

@testset "Spectral Spreading - JONSWAP" begin
    # 1. Create a JONSWAP continuous spectrum
    Hs = 3.0
    Tp = 10.0
    m0_exact = (Hs / 4.0)^2
    spec = JONSWAP(Hs, Tp)
    
    # 2. Test: Check if Ag was calculated
    @test spec.Ag > 0
    
    # 3. Test: Integration of the density
    # The area under S(f) should be (Hs/4)^2 after normalization by Ag
    # m0 = ∫ S(f) df
    m0_num = gaussQuad1D_4(f -> get_density(spec, f), get_fmin(spec), get_fmax(spec), 500)
    Hs_calc = 4 * sqrt(m0_num)
    
    @info "Target Hs: $Hs, Calculated Hs: $Hs_calc"
    @test Hs_calc ≈ Hs rtol = 1e-2
    @info "Target m0: $m0_exact, Numerical m0:   $m0_num"
    @test m0_num ≈ m0_exact rtol = 1e-2

    # 4. Test: Check that fmax is sufficient to integrate 99% of area under the curve (contains almost all tail)
    m0_ext = gaussQuad1D_4(f -> get_density(spec, f), get_fmin(spec), get_fmax(spec, multiplier = 20.0), 1000)
    energy_ratio = m0_num / m0_ext

    @info "Default max frequency: $get_fmax(spec)"
    @info "Energy Coverage: $(round(energy_ratio * 100, digits=4))%"
    @test energy_ratio >= 0.99
end

end # module