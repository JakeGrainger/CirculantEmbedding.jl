function approx_cov(Γ::Kernel{D,P}, lags) where {D,P}
    nfreq = choose_nfreq.(lags.iterators)
    Δ = step.(lags.iterators)
    n = length.(lags.iterators)
    freq = fftfreq.(nfreq, 1.0./Δ)

    f = aliased_sdf.(Γ, Ref(Iterators.product(freq...)), Ref(Δ))
    c = FFTView(fft_array(f))
    return real.(c[CartesianIndices(ntuple(i->-n[i]:n[i]-1,Val{D}()))]) .* prod(step, freq)
end

choose_nfreq(lag) = max(length(lag), 2048)
