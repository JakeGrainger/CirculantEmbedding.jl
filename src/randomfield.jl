abstract type RandomField{D,P} end

struct IndependentFields{D,P,T<:NTuple{P,<:RandomField{D,1}}} <: RandomField{D,P}
    fields::T
    function IndependentFields(fields::T) where {D,P,T<:NTuple{P,<:RandomField{D,1}}}
        @assert all(x->getmesh(x)==getmesh(fields[1]), fields)
        new{D,P,T}(fields)
    end
end
getmesh(f::IndependentFields) = getmesh(f.fields[1])

function rand(f::IndependentFields{D,P,T}) where {D,P,T}
    X = rand.(f.fields)
    return [SVector{P,Float64}(X[p][ind] for p in 1:P) for ind in CartesianIndices(X[1])]
end

Distributions.mean(f::IndependentFields{D,P,T}) where {D,P,T} = SVector(getindex.(mean.(f.fields),1))
Distributions.cov(f::IndependentFields{D,P,T}, lag) where {D,P,T} = SMatrix{P,P,Float64,P*P}(i==j ? cov(f.fields[i],lag) : 0.0 for i in 1:P, j in 1:P)
Distributions.var(f::IndependentFields{D,P,T}) where {D,P,T} = SMatrix{P,P,Float64,P*P}(i==j ? var(f.fields[i]) : 0.0 for i in 1:P, j in 1:P)
sdf(f::IndependentFields{D,P,T},freq) where {D,P,T} = SMatrix{P,P,ComplexF64,P*P}(i==j ? sdf(f.fields[i],freq) : 0.0 for i in 1:P, j in 1:P)