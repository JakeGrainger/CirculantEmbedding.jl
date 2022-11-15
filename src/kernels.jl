abstract type Kernel{D,P} end

Base.broadcastable(k::Kernel) = Ref(k)
Distributions.var(k::Kernel{D,P}) = cov(k, @SVector zeros(D))

struct AdditiveKernel{D,P,K1<:Kernel{D,P},K2<:Kernel{D,P}} <: Kernel{D,P}
    k1::K1
    k2::K2
end
Base.:+(k1::Kernel,k2::Kernel) = AdditiveKernel(k1,k2)
Distributions.cov(k::AdditiveKernel, h) = cov(k.k1, h) + cov(k.k2, h)
sdf(k::AdditiveKernel, w) = sdf(k.k1, w) + sdf(k.k2, w)

include("kernels/Matern.jl")
include("kernels/OscillatoryMatern.jl")