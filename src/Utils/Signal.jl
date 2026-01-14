module Signal 

"""
FFT
===========
Module providing tools for spectral reconstruction through signal theory methodology.
"""

using FFTW
using ..AngularSpreading
using ..ContinuousSpectrums
using ..FrequencySampling
using ..SpectralSpreading
using ..AiryWaves

export generate_signal, get_single_sided_spectrum, averaged_psd


function generate_signal(state::AiryState, N::Int; fs::Real = 2.0 *get_fmax(state.spectrum.spectrum), P::AbstractVector{<:Real} = [0.0, 0.0, 0.0], vars = [:η]) 
    """
    Generates signal (time serie) corresponding at airy state "state" and wave profile specified by "vars" at point "P" and constituted of "N" samples with frequency "fs".

    Inputs:
      - state : airy waves structure containing discrete spectrum and angular spreading model
      - N     : number of samples (should be an even number, power of 2 for instance)
      - fs    : window sampling frequency (default = 2*fmax -> Nyquist frequency -> sampling frequency fs must be MORE than double the signal maximum frequency -> AVOID ALIASING)
      - P     : spatial coordinates in airy sea where signal is generated
      - vars  : airy profile from which generate the signal (default is the sea elevation)
    """

    # Time array
    dt = 1.0 / fs
    t = collect(0:dt:(N-1)*dt)

    # Generate signal
    profiles = sea_profiles(state, [P[1]], [P[2]], [P[3]], t, vars=vars)

    # Extract the signal
    target_var = vars isa Symbol ? vars : vars[1]
    return t, profiles[target_var][1,1,1,:]
end


function generate_signal(discspectrum::DiscreteSpectrum, h::Real, N::Int; fs::Real = 2.0 *get_fmax(discspectrum.spectrum), P::AbstractVector{<:Real} = [0.0, 0.0, 0.0], vars = [:η])
    """
    Generates signal (time serie) corresponding at airy state with discrete spectrum "discspectrum" and wave profile specified by "vars" at point "P" and constituted of "N" samples with frequency "fs".

    Inputs:
      - discspectrum : airy waves structure containing discrete spectrum and angular spreading model
      - h            : water depth (m)
      - N            : number of samples (should be an even number, power of 2 for instance)
      - fs           : window sampling frequency (default = 2*fmax -> Nyquist frequency -> sampling frequency fs must be MORE than double the signal maximum frequency -> AVOID ALIASING)
      - P            : spatial coordinates in airy sea where signal is generated
      - vars         : airy profile from which generate the signal (default is the sea elevation)
    """
    
    # Generate null angular spreading structure
    AngSprea = SpreadingModel(0.0)

    # Generate Airy sea state
    airysea = AiryState(discspectrum, AngSprea, h)

    # Generate signal
    return generate_signal(airysea, N, fs = fs, P = P, vars = vars)
end


function get_single_sided_spectrum(signal::AbstractVector{<:Real}, fs::Real)
    """
    Compute the single-sided amplitude and power spectrum of a "signal" with sampling frequency "fs".
 
    Input
    ----------
    - signal : Input signal samples (1D array).
    - fs     : Sampling frequency of the input signal in Hz.
 
    Returns
    -------
    - fHalf  : Array of frequency (Hz) for the single-sided spectrum.
    - fAmp   : Single-sided amplitude spectrum of the input signal.
    - psd    : Single-sided power spectral density (PSD) of the input signal.
 
    Notes
    -----
    - The function ensures the input signal length is even for FFT computation.
    - The amplitude spectrum is normalized and scaled for single-sided representation.
    - The power spectral density is computed per frequency bin.
    - Also prints sample length, frequency resolution, and maximum frequency.
 
    """
 
    N = length(signal) - mod(length(signal), 2)   # even size of signal samples
    signal = signal[1:N]                          # truncate signal
 
    println("Sample Length = ", N)
    println("Least count Hz = ", fs / N)
    println("Max Freq (Half band) Hz = ", fs / 2)
 
    # frequencies of half spectrum
    fHalf = fs * range(0, N/2 - 1) / N

    # Raw FFT magnitude
    fAmp = fft(signal)
    fAmp = 2.0 * abs.(fAmp[1:Int(N/2)]) / N

    # Power spectral density
    psd = fAmp .^ 2 / 2.0 / (fHalf[2] - fHalf[1])
 
    return fHalf, fAmp, psd
end


function averaged_psd(fHalf, psd, navg=20)
    """
    Averages power spectral density over multiple samples.
    """
    n = length(psd)
    n_smooth = div(n, navg)
    f_res = zeros(n_smooth)
    psd_res = zeros(n_smooth)
    
    for i in 1:n_smooth
        idx = ((i-1)*navg + 1):(i*navg)
        f_res[i] = mean(fHalf[idx])
        psd_res[i] = mean(psd[idx])
    end
    return f_res, psd_res
end

end