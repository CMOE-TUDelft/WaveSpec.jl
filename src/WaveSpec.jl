module WaveSpec

# 1. Include the foundational utilities first
include("Utils/Truncation.jl")
@reexport using .Truncation

# 2. Include the Spreading modules (which depend on Truncation)
include("SpectralSpreading.jl")
@reexport using .SpectralSpreading

include("AngularSpreading/AngularSpreading.jl")
@reexport using .AngularSpreading

# 3. Include the Physics/Synthesis module (which depends on the above)
include("AiryWaves.jl")
@reexport using .AiryWaves


"""
    WaveSpec
A Julia package for generating stochastic sea states using Airy Wave theory,
customizable spectral shapes (JONSWAP, Pierson-Moskowitz), and angular spreading models.
"""
function __init__()
    @info "WaveSpec initialized: Ready for spectral sea state synthesis."
end

end