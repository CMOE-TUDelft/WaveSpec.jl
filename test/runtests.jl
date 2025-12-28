module WaveSpecTests

using Test

@time @testset "Continuous Spectrums" begin include("ContinuousSpectrumsTests/runtests.jl") end

@time @testset "Discrete Spectrums" begin include("DiscreteSpectrumsTests/runtests.jl") end

@time @testset "Angular Spreading Models" begin include("AngularSpreadingTests/runtests.jl") end

@time @testset "Truncation" begin include("TruncationTests.jl") end

end # module