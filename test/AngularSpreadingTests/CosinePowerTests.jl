module CosinePowerTests

using Test
using WaveSpec.AngularSpreading
using WaveSpec.Truncation
using Distributions

# Define a Cosine Power Spreading Model for testing
μ = 0.0
n = 4
sampling = CosinePowerDistribution(μ, n)
# Truncate to a full circle
model = TruncatedModel(sampling, -3*π/4, 3*π/4)

@testset "Truncated Cosine Power Spreading Specifics Tests" begin
    # 1. Check physical bounds of Cosine Power
    # Density should be 0 at ±π/2 even if truncation is at ±π
    @test pdf(model, 4*π/5) ≈ 0.0 atol=1e-7
    @test pdf(model, -4*π/5) ≈ 0.0 atol=1e-7
    
    # 2. Check Peak
    @test pdf(model, μ) > 0.0
end

end # module CosinePowerTests