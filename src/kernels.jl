abstract type Kernel{D,P} end

Base.broadcastable(k::Kernel) = Ref(k)

struct Matern{D,P,T,L} <: Kernel{D,P}
    σ²::SHermitianCompact{P,T,L}
    ν::SHermitianCompact{P,T,L}
    a::SHermitianCompact{P,T,L}
    function Matern(σ²::SHermitianCompact{P,T,L},ν::SHermitianCompact{P,T,L},a::SHermitianCompact{P,T,L},::Val{D}) where {D,P,T,L}
        new{D,P,T,L}(σ²,ν,a)
    end
    Matern(σ²,ν,a,D=1) = Matern(σ²,ν,a,Val{D}())
end


function Distributions.cov(Γ::Matern{D,1,T,L}, h) where {D,T,L}
    nh = norm(h)
    Γ.σ²[1,1] * materncorr(Γ.ν[1,1], Γ.a[1,1], nh)
end

function Distributions.cov(Γ::Matern{D,P,T,L}, h) where {D,P,T,L}
    nh = norm(h)
    return SMatrix{P,P,T,P^2}(Γ.σ²[i,j] * materncorr(Γ.ν[i,j], Γ.a[i,j], nh) for i in 1:P, j in 1:P)
end

function materncorr(ν, a, nh)
    nh < 1e-10 ? 1.0 : 2^(1-ν)/gamma(ν) * (a*nh)^ν * besselk(ν, a*nh)
end

function sdf(Γ::Matern{D,P,T,L}, w) where {D,P,T,L}
    nw = norm(w)
    return SMatrix{P,P,complex(T),P^2}(Γ.σ²[i,j] * materncorr_ft(Γ.ν[i,j], Γ.a[i,j], nw, D) for i in 1:P, j in 1:P)
end

function materncorr_ft(ν, a, nw, d)
    2^d * π^(d/2)   * gamma(ν+d/2)/gamma(ν) * a^(2ν) * (1/(a^2+4*π^2*nw^2))^(ν+d/2)
end