module CirculantEmbedding

    using FFTW, StaticArrays, Meshes, SpecialFunctions, Distributions, LinearAlgebra, FFTViews
    import Base: rand

    include("kernels.jl")
    include("simulation.jl")
    include("fft_array.jl")
    include("cov_approximation.jl")

    export GaussianProcess, Matern, OscillatoryMatern, sdf
end
