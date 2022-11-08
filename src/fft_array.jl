function fft_array(args...;kwargs...) # default fallback method when scalar
    fft(args..., kwargs...)
end
function fft_array!(args...;kwargs...) # default fallback method when scalar
    fft!(args..., kwargs...)
end
 
function fft_array(x::Array{A,N},dims=1:ndims(x);kwargs...) where {A<:AbstractArray,N}
    x_out = alloc_fft_output(x)
    fft_array!(x_out, dims;kwargs...)
end
function fft_array!(x::Array{A,N},dims=1:ndims(x);kwargs...) where {A<:AbstractArray,N}
    checksamesize(x)
    xflat = reinterpret(reshape, deep_eltype(x), x)
    fft!(xflat,1 .+ dims;kwargs...)
    return x
end

function ifft_array(args...;kwargs...) # default fallback method when scalar
    ifft(args..., kwargs...)
end
function ifft_array!(args...;kwargs...) # default fallback method when scalar
    ifft!(args..., kwargs...)
end
 
function ifft_array(x::Array{A,N},dims=1:ndims(x);kwargs...) where {A<:AbstractArray,N}
    x_out = alloc_fft_output(x)
    ifft_array!(x_out, dims;kwargs...)
end
function ifft_array!(x::Array{A,N},dims=1:ndims(x);kwargs...) where {A<:AbstractArray,N}
    checksamesize(x)
    xflat = reinterpret(reshape, deep_eltype(x), x)
    ifft!(xflat,1 .+ dims;kwargs...)
    return x
end
 
deep_eltype(x) = deep_eltype(typeof(x))
deep_eltype(x::Type{<:AbstractArray}) = deep_eltype(eltype(x))
deep_eltype(x::Type) = x
 
function checksamesize(x::Array) # fallback
    all(size(x[i])==size(x[1]) for i in eachindex(x)) || throw(DimensionMismatch("All arrays in array should be the same type."))
    return nothing
end
function checksamesize(::Array{T,N}) where {T<:StaticArray{S,T₂,D},N} where {S,T₂,D}
    return nothing # no check necessary in this case.
end

function alloc_fft_output(x::Array{A,D}) where {A,D}
    y=Array{mycomplex(A),D}(undef, size(x))
    for i in eachindex(y,x)
        y[i] = mycomplex(x[i])
    end
    return y
end

mycomplex(x) = Base.complex(x)
mycomplex(::Type{SArray{NTuple{D,S},T,D,N}}) where {D,T,S,N} = SArray{NTuple{D,S},mycomplex(T),D,N}
mycomplex(::Type{SHermitianCompact{D, T, L}}) where {D,T,L} = SHermitianCompact{D, mycomplex(T), L}
mycomplex(s::SHermitianCompact{D, T, L}) where {D, T<:Real, L} = SHermitianCompact(mycomplex(s.lowertriangle))