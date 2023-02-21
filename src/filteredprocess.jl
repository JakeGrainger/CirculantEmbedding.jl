struct FilteredRandomField{D,P,Q,L<:RandomField{D,Q},T<:LinearFilter{Q,P,Float64}} <: RandomField{D,P}
    latent::L
    transform::T
end

function rand(f::FilteredRandomField, show_latent=true)
    latent_fields = rand(f.latent)
    output_fields = f.transform.(latent_fields)
    if show_latent
        return (output_fields, latent_fields)
    else
        return output_fields
    end
end