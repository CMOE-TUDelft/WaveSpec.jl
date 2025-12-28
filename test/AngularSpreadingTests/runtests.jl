module AngularSpreadingTests

using Test

@testset "Cosine Power Model" begin include("CosinePowerTests.jl") end

@testset "Von Mises Model" begin include("VonMisesTests.jl") end

@testset "Normalization" begin include("NormalizationTests.jl") end

end # module