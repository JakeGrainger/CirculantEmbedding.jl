struct IsoOscillatoryMatern{D,P,T,L} <: KernelSdfOnly{D,P}
    matern::Matern{D,P,T,L}
    offset::T
end
IsoOscillatoryMatern(σ²,ν,a,offset,D=1) = IsoOscillatoryMatern(Matern(σ²,ν,a,Val{D}()), offset)

function sdf(k::IsoOscillatoryMatern, h)
    return sdf(k.matern, abs(norm(h)-k.offset))
end