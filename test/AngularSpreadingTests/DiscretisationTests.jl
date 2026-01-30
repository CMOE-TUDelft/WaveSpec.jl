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
    models = [("Normal", DiscreteAngularSpreading(:normal, μ, σ, θmin, θmax, nθ)),
              ("Uniform", DiscreteAngularSpreading(:uniform, μ, σ, θmin, θmax, nθ)),
              ("Cosine Power", DiscreteAngularSpreading(:cosinepow, μ, σ, θmin, θmax, nθ)),
              ("Von Mises", DiscreteAngularSpreading(:vonmises, μ, σ, θmin, θmax, nθ)),
              ("Donelan-Banner", DiscreteAngularSpreading(:donelan, μ, σ, θmin, θmax, nθ))]
    return models
end


function test_model_discretisation(model::DiscreteAngularSpreading; rtol=1e-3)
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

@testset "get_bandwidth bounds checking" begin
    # Create a simple test model
    model = DiscreteAngularSpreading(:normal, 0.0, 15.0 * (π/180), -π/3, π/3, 37)
    
    # Test valid index returns Float64
    @test typeof(AngularSpreading.get_bandwidth(model, 1)) == Float64
    @test typeof(AngularSpreading.get_bandwidth(model, model.nθ - 1)) == Float64
    
    # Test invalid indices throw BoundsError
    @test_throws BoundsError AngularSpreading.get_bandwidth(model, 0)
    @test_throws BoundsError AngularSpreading.get_bandwidth(model, -1)
    @test_throws BoundsError AngularSpreading.get_bandwidth(model, model.nθ)
    @test_throws BoundsError AngularSpreading.get_bandwidth(model, model.nθ + 100)
end

end # module DiscretisationTests