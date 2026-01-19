module AngularSpreadingTests

using Test

@testset "Continuous Angular Spreading Specific Tests" begin 
    include("CosinePowerTests.jl") 
    include("VonMisesTests.jl")
end

@testset "Discrete Angular Spreading Models Tests" begin 
    include("TruncationTests.jl")
    include("DiscretisationTests.jl") 
end

end # module