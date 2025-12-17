module FrequencySampling

export AbstractSampling, UniformSampling, LogSampling, SpectralSampling

abstract type AbstractSampling end

struct UniformSampling <: AbstractSampling end
struct LogSampling <: AbstractSampling end
struct SpectralSampling <: AbstractSampling end

"""
    generate_grid(strategy, spectrum, f_min, f_max, n_f)

Dispatches the grid generation based on the selected strategy.
"""
function generate_grid(::UniformSampling, _, fmin, fmax, n_f)
    return collect(range(fmin, fmax, length=n_f))
end

function generate_grid(::LogSampling, _, fmin, fmax, n_f)
    return exp.(range(log(fmin), log(fmax), length=n_f))
end


function generate_grid(::SpectralSampling, spectrum_type, fmin, fmax, n_f)
    
end

end