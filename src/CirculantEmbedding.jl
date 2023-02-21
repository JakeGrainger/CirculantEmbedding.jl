module CirculantEmbedding

    using FFTW, StaticArrays, Meshes, SpecialFunctions, Distributions, LinearAlgebra, FFTViews
    import Base: rand

    include("kernels.jl")
    include("randomfield.jl")
    include("gaussianprocess.jl")
    include("linearfilter.jl")
    include("filteredprocess.jl")
    include("fft_array.jl")
    include("cov_approximation.jl")

    export GaussianProcess, Matern, OscillatoryMatern, sdf, approx_cov, IndependentFields, LagZeroFilter, FilteredRandomField, getmesh
end
