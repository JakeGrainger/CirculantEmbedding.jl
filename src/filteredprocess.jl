struct FilteredRandomField{D,P,Q,L<:RandomField{D,Q},T<:LinearFilter{Q,P,Float64}} <: RandomField{D,P}
    latent::L
    transform::T
end

function rand(f::FilteredRandomField)
    latent_fields = rand(f.latent)
    output_fields = f.transform(latent_fields)
    return (rf = output_fields, latent = latent_fields)
end

getmesh(f::FilteredRandomField) = getmesh(f.latent)

Distributions.mean(f::FilteredRandomField{D,P,Q,L,LagZeroFilter{Q,P,Float64,S}}) where {D,P,Q,L,S} = f.transform.A * mean(f.latent)
Distributions.cov(f::FilteredRandomField{D,P,Q,L,LagZeroFilter{Q,P,Float64,S}}, lag) where {D,P,Q,L,S} = f.transform.A * cov(f.latent,lag) * f.transform.A'
Distributions.var(f::FilteredRandomField{D,P,Q,L,LagZeroFilter{Q,P,Float64,S}}) where {D,P,Q,L,S} = f.transform.A * var(f.latent) * f.transform.A'
sdf(f::FilteredRandomField{D,P,Q,L,LagZeroFilter{Q,P,Float64,S}}, freq) where {D,P,Q,L,S} = f.transform.A * sdf(f.latent, freq) * f.transform.A'