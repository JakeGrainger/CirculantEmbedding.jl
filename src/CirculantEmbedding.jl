module CirculantEmbedding

    using FFTW, StaticArrays, Meshes, SpecialFunctions, Distributions, LinearAlgebra
    include("kernels.jl")
    include("simulation.jl")
    include("fft_array.jl")

    export GaussianProcess, Matern, simulate
end
