struct FilteredRandomField{D,P,Q,L<:RandomField{D,Q},T<:LinearFilter{Q,P,Float64}} <: RandomField{D,P}
    latent::L
    transform::T
end

function rand(f::FilteredRandomField)
    latent_fields = rand(f.latent)
    output_fields = f.transform.(latent_fields)
    return (rf = output_fields, latent = latent_fields)
end

getmesh(f::FilteredRandomField) = getmesh(f.latent)