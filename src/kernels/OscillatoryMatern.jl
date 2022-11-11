struct OscillatoryMatern{D,P,T,L} <: Kernel{D,P}
    matern::Matern{D,P,T,L}
    offset::SVector{D,T}
end
OscillatoryMatern(σ²,ν,a,offset::SVector{D,T}) where {D,T} = OscillatoryMatern(Matern(σ²,ν,a,Val{D}()), offset)
OscillatoryMatern(σ²,ν,a,offset::NTuple{D,T}) where {D,T}  = OscillatoryMatern(σ²,ν,a,SVector{D,T}(offset))

Distributions.cov(k::OscillatoryMatern, h::NTuple{D,T}) where {D,T} = cov(k, SVector{D,T}(h))
Distributions.cov(k::OscillatoryMatern, h::SVector) = cos(2π*(k.offset'*h))*cov(k.matern,h)

sdf(k::OscillatoryMatern, h::NTuple{D,T}) where {D,T} = sdf(k, SVector{D,T}(h))
function sdf(k::OscillatoryMatern, h::SVector)
    return (sdf(k.matern, h-k.offset)+sdf(k.matern, h+k.offset))/2
end