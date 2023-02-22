"""
    LinearFilter{Q,P,T}

A linear filter from a Q dimensional space to a P dimensional space, with eltype T.
"""
abstract type LinearFilter{Q,P,T} end

struct LagZeroFilter{Q,P,T,L} <: LinearFilter{Q,P,T}
    A::SMatrix{P,Q,T,L}
end

function (f::LagZeroFilter{Q,P,T,L})(X::Array{SVector{Q,S},D}) where {D,Q,P,T,L,S}
    return map(x->f.A*x, X)
end

transferfunction(f::LagZeroFilter, freq) = f.A