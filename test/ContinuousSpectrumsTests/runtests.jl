module ContinuousSpectrumsTests

using Test

# --- Specific Spectrum Tests ---

@testset "Continuous Spectrum Specific Tests" begin 
    include("JONSWAPTests.jl") 
    include("TMATests.jl")
    include("DonelanTests.jl")
    include("BretschneiderTests.jl")
    include("OchiHubbleTests.jl")
end

# --- Integration Tests for all Spectra ---

@testset "Spectral Integration Tests" begin include("IntegrationTests.jl") end

end # module