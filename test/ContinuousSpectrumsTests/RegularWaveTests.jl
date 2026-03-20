using Test
using WaveSpec.ContinuousSpectrums
import WaveSpec.SpectralSpreading as Spreading
import WaveSpec.SpectralSampling as Sampling

@testset "RegularWave" begin
    H = 2.0
    T = 8.0
    s = RegularWave(H, T)
    @test get_Hs(s) ≈ H
    @test get_Tp(s) ≈ T
    f0 = 1 / T
    S0 = get_density(s, f0)
    # Should be much larger at f0 than away from f0
    @test S0 > 0.1
    @test get_density(s, f0 + 0.1) < S0
    ds = Spreading.DiscreteSpectralSpreading(s, Sampling.UniformSampling(),0.01, 0.5, 100)
    area = sum(Spreading.get_integrated_energies(ds))
    @test area ≈ H^2/8 atol=1e-2
end
