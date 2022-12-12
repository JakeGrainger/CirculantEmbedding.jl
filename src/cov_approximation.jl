"""
    approx_cov(Γ::Kernel{D,P}, lags) where {D,P}

Function to approximate the covariance from the sdf.
"""
function approx_cov(Γ::Kernel{D,P}, lags) where {D,P}
    nfreq = choose_nfreq.(lags.iterators)
    Δ = step.(lags.iterators)
    n = length.(lags.iterators) .÷ 2
    freq = fftfreq.(nfreq, 1.0./Δ)

    f = aliased_sdf.(Γ, Iterators.product(freq...), Ref(Δ))
    c = FFTView(fft_array(f))
    return fftshift(real.(c[CartesianIndices(ntuple(i->-n[i]:n[i]-1,Val{D}()))]) .* prod(step, freq))
end

"""
    choose_nfreq(lag)

Chooses the maximum number of frequencies to use to approximate a given number of lags.
"""
choose_nfreq(lag) = max(length(lag), 2048)
