abstract type RandomField{D,P} end

struct IndependentFields{D,P,T<:NTuple{P,<:RandomField{D,1}}} <: RandomField{D,P}
    fields::T
end

function rand(f::IndependentFields{D,P,T}) where {D,P,T}
    X = rand.(f.fields)
    return [SVector{P,Float64}(X[p][ind] for p in 1:P) for ind in eachindex(X...)]
end

Distributions.mean(f::IndependentFields{D,P}) where {D,P} = mean.(f.fields)
Distributions.cov(f::IndependentFields{D,P}, lag) where {D,P} = SMatrix{P,P,Float64,P*P}(i==j ? cov(f.fields[i],lag) : 0.0 for i in 1:P, j in 1:P)
Distributions.var(f::IndependentFields{D,P}) where {D,P} = SMatrix{P,P,Float64,P*P}(i==j ? var(f.fields[i]) : 0.0 for i in 1:P, j in 1:P)
sdf(f::IndependentFields{D,P},freq) where {D,P} = SMatrix{P,P,ComplexF16,P*P}(i==j ? sdf(f.fields[i],freq) : 0.0 for i in 1:P, j in 1:P)