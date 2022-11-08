module CirculantEmbedding

    using FFTW, StaticArrays, Meshes, SpecialFunctions, Distributions, LinearAlgebra
    import Base: rand

    include("kernels.jl")
    include("simulation.jl")
    include("fft_array.jl")

    export GaussianProcess, Matern, simulate
end
