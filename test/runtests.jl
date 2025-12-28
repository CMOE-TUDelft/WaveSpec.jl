module WaveSpecTests

using Test

@time @testset "Continuous Spectrums" begin include("ContinuousSpectrums/runtests.jl") end

@time @testset "Discrete Spectrums" begin include("test_DiscreteSpectrum.jl") end

@time @testset "Angular Spreading Models" begin include("AngularSpreadingTests/runtests.jl") end

@time @testset "Truncation" begin include("TruncationTests.jl") end

end # module