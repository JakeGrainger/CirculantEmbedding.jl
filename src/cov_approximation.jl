function approx_cov(Γ::Kernel{D,P}, lags) where {D,P}
    nfreq = choose_nfreq.(lags)
    Δ = step.(lags)
    n = length.(lags)
    freq_f = fftfreq.(nfreq, 1.0./Δ)

    f = [aliased_sdf(Γ, (l1,l2)) for l1 in freq_f, l2 in freq_f]
    c = FFTView(fft_array(f))
    return real.(c[CartesianIndices(ntuple(i->-n[i]:n[i]-1,Val{D}()))]) .* (step(freq_f))^2
end

choose_nfreq(lag) = max(length(lag), 2048)
