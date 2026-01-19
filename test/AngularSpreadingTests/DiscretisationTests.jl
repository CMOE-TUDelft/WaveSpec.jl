module DiscretisationTests

using Test
using Distributions
using WaveSpec.AngularSpreading
using WaveSpec.Integration

"""
    prepare_discrete_models()
Prepares a list of implemented discrete truncated angular models, with tuples: (model_name::String, model_instance::AbstractSpectrum)
"""
function prepare_discrete_models()
    μ = 0.0
    σ = 15.0 * (π/180)
    θmin = -π/3
    θmax = π/3
    nθ = 37

    # Prepare Models in list of tuples (name, structure)
    models = [("Normal", SpreadingModel(:normal, μ, σ, θmin, θmax, nθ)),
              ("Uniform", SpreadingModel(:uniform, μ, σ, θmin, θmax, nθ)),
              ("Cosine Power", SpreadingModel(:cosinepow, μ, σ, θmin, θmax, nθ)),
              ("Von Mises", SpreadingModel(:vonmises, μ, σ, θmin, θmax, nθ)),
              ("Donelan-Banner", SpreadingModel(:donelan, μ, σ, θmin, θmax, nθ))]
    return models
end


function test_model_discretisation(model::SpreadingModel; rtol=1e-3)
    # Discrete model area must be 1.0 
    weights = get_weights(model)
    Δθ = get_bandwidths(model)
    area = sum(weights .* Δθ)
    
    @test area ≈ 1.0 rtol=rtol
end


@testset "Angular Spreading Discretisation Tests" begin
    models = prepare_discrete_models()
    for (name, model) in models
        @testset "Model: $name" begin
            test_model_discretisation(model)
        end
    end
end

end # module DiscretisationTests