abstract type Kernel{D,P} end

Base.broadcastable(k::Kernel) = Ref(k)
"""
    var(k::Kernel{D,P}) where {D,P}

Compute the variance of a stationary kernel.
"""
Distributions.var(k::Kernel{D,P}) where {D,P} = cov(k, @SVector zeros(D))

struct AdditiveKernel{D,P,K1<:Kernel{D,P},K2<:Kernel{D,P}} <: Kernel{D,P}
    k1::K1
    k2::K2
end
Base.:+(k1::Kernel,k2::Kernel) = AdditiveKernel(k1,k2)
Distributions.cov(k::AdditiveKernel, h) = cov(k.k1, h) + cov(k.k2, h)
sdf(k::AdditiveKernel, w) = sdf(k.k1, w) + sdf(k.k2, w)

abstract type KernelSdfOnly{D,P} <: Kernel{D,P} end
Distributions.cov(k::KernelSdfOnly, h) = error("Covariance not defined for $(typeof(k)). Can still be approximated, see approx_cov.")

"""
    aliased_sdf(Γ::Kernel{D,P}, k, Δ; K = 3) where {D,P}

Compute the aliased sdf for a given kernel.
"""
function aliased_sdf(Γ::Kernel{D,P}, k, Δ; K = 3) where {D,P}
    return sum(sdf(Γ, k .+ j.I ./ Δ) for j in CartesianIndices(ntuple(j->-K:K, Val{D}())))
end

include("kernels/Matern.jl")
include("kernels/OscillatoryMatern.jl")
include("kernels/IsotropicOscillatoryMatern.jl")