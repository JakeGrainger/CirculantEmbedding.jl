abstract type Kernel{D,P} end

Base.broadcastable(k::Kernel) = Ref(k)

include("kernels/Matern.jl")
include("kernels/OscillatoryMatern.jl")