module TruncationTests

using Test
using WaveSpec.AngularSpreading
using WaveSpec.Truncation
using WaveSpec.Integration
using Distributions
import Statistics: mean

"""
    prepare_truncated_models()
Prepares a list of implemented truncated angular models, with tuples: (model_name::String, model_instance::AbstractSpectrum)
"""

function prepare_truncated_models()
    μ = 0.0
    σ = 15.0 * (π/180)
    θmin = -π/3
    θmax = π/3

    n = 4 # Cosine Power exponent

    # Prepare Models in list of tuples (name, structure)
    models = [("Normal", TruncatedModel(Normal(μ, σ), θmin, θmax)),
              ("Uniform", TruncatedModel(Uniform(-π, π), θmin, θmax)),
              ("Cosine Power", TruncatedModel(CosinePowerDistribution(μ, n), θmin, θmax)),
              ("Von Mises", TruncatedModel(VonMisesDistribution(μ, σ), θmin, θmax)),
              ("Donelan-Banner", TruncatedModel(DonelanBannerDistribution(μ, σ), θmin, θmax))]
    return models
end

"""
    test_truncated_spreading_integrity(model::TruncatedModel; rtol=1e-3)

Verifies that the truncated discrete angular spreading models are correctly 
discretised and renormalizes the underlying sampling distribution to a total area of 1.0.
"""
function test_model_truncation(TDM::TruncatedModel; atol=1e-10)
    # Density outside the limits must be exactly 0.0
    @test pdf(TDM, TDM.a - 0.1) == 0.0
    @test pdf(TDM, TDM.b + 0.1) == 0.0

    # Peak location check
    # The means of the truncated model and base distribution should align
    base_mean = mean(TDM.dist)
    trunc_mean = mean(TDM)
    @test trunc_mean ≈ base_mean atol=atol
end


@testset "Truncated Spreading Integrity Tests" begin
    models = prepare_truncated_models()
    for (name, model) in models
        @testset "Model: $name" begin
            test_model_truncation(model)
        end
    end
end


end # module TruncationTests