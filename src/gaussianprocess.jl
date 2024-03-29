struct GaussianProcess{M<:Real,D,P,T<:Real,G<:Mesh{D,T},S} <: RandomField{D,P}
    mean::SVector{P,M}
    Γ::Kernel{D,P}
    mesh::G
    sim_prealloc::S
end
"""
    GaussianProcess(mean, Γ::Kernel, mesh::Mesh; pad=0)

Construct a Gaussian process with kernel `Γ` on a mesh `mesh`.
"""
GaussianProcess(μ::SVector{P,M}, Γ::Kernel{D,P}, mesh::Mesh; pad=0) where {D,P,M} = GaussianProcess(μ, Γ, mesh, preallocate_gp_simulation(Γ, mesh, pad))
GaussianProcess(μ::M, Γ::Kernel{D,P}, mesh::Mesh; pad=0) where {D,P,M<:Real} = GaussianProcess(SVector{P,M}(ntuple(d->μ, Val{P}())), Γ, mesh; pad=pad)
GaussianProcess(Γ::Kernel, mesh::Mesh; pad=0) = GaussianProcess(0.0, Γ, mesh; pad=pad)

getmesh(g::GaussianProcess) = g.mesh
Distributions.mean(g::GaussianProcess) = g.mean
Distributions.var(g::GaussianProcess) = Distributions.var(g.Γ)
Distributions.cov(g::GaussianProcess, h) = Distributions.cov(g.Γ, h)
sdf(g::GaussianProcess, freq) = sdf(g.Γ, freq)

preallocate_gp_simulation(Γ::Kernel, mesh::Mesh, pad) = error("Simulation only currently implemented on a regular grid.")
preallocate_gp_simulation(Γ::Kernel{D,P}, mesh::CartesianGrid{D,T}, pad) where {D,P,T} = CirculantPrealloc(Γ, mesh, pad)
preallocate_gp_simulation(Γ::Kernel{D,P}, mesh::CartesianGrid{D2,T}, pad) where {D,D2,P,T} = error("Kernel has different dimension to grid.")

struct CirculantPrealloc{T,S,D}
    L::Array{T,D}
    Y::Array{S,D}
    function CirculantPrealloc(Γ::Kernel{D,P}, mesh::CartesianGrid{D,T}, pad=0) where {D,P,T}
        C = SHermitianCompact.(cov_at_ft(Γ, mesh, pad)) # strictly not symmetric, but doesn't change the algorithm as result of fft will be hermitian symm.
        A = fft_array(C)
        L = cholesky.(A)
        Y = Array{SVector{P,ComplexF64},D}(undef, size(L))
        all(x->!isnan(sum(x)), L) || error("Some NaNs found, check parameters.")
        new{eltype(L), eltype(Y), D}(L,Y)
    end
    function CirculantPrealloc(Γ::Kernel{D,1}, mesh::CartesianGrid{D,T}, pad=0) where {D,T}
        C = cov_at_ft(Γ, mesh, pad)
        A = fft(C)
        L = sqrt.(A)
        Y = Array{ComplexF64,D}(undef, size(L))
        all(!isnan, L) || error("Some NaNs found, check parameters.")
        new{eltype(L), eltype(Y), D}(L,Y)
    end
end

function cov_at_ft(Γ::Kernel, mesh, pad=0)
    lags = compute_lags(mesh, pad)
    return multi_cov(Γ, lags)
end

const KernelSomeSdfOnly = Union{KernelSdfOnly, AdditiveKernel{<:KernelSdfOnly,<:Kernel}, AdditiveKernel{<:Kernel, <:KernelSdfOnly}, AdditiveKernel{<:KernelSdfOnly, <:KernelSdfOnly}}
function multi_cov(Γ::Kernel, lags)
    return cov.(Γ, lags)
end
function multi_cov(Γ::Union{KernelSdfOnly, AdditiveKernel}, lags)
    return approx_cov(Γ, lags)
end

"""
    compute_lags(mesh::CartesianGrid{D,T}, pad) where {D,T}

Compute lags required for circulant embedding on a given mesh.
"""
function compute_lags(mesh::CartesianGrid{D,T}, pad) where {D,T}
    n, Δ = mesh.topology.dims, mesh.spacing
    m = ntuple(d -> choose_circ_m(n[d], pad), Val{D}())
    Iterators.product(ntuple(d -> fftfreq(m[d],m[d]*Δ[d]), Val{D}())...)
end

"""
    choose_circ_m(n::Int, pad)

Choose the `m` used in circulant embedding.
"""
function choose_circ_m(n::Int, pad)
    return 2n-2+pad
end

function rand(gp::GaussianProcess{M,D,P,T,CartesianGrid{D,T},S}) where {M,D,P,T,S}
    for i in eachindex(gp.sim_prealloc.L,gp.sim_prealloc.Y)
        Z = SVector{P,complex(T)}(complex(rand(Normal(0,1)), rand(Normal(0,1))) for i in 1:P)
        gp.sim_prealloc.Y[i] = gp.sim_prealloc.L[i].L * Z
    end
    W = fft_array(gp.sim_prealloc.Y) / sqrt(length(gp.sim_prealloc.Y))
    return gp.mean .+ real.(W[CartesianIndices(size(gp.mesh))])#, imag.(W[CartesianIndices(size(gp.mesh))])
end
function rand(gp::GaussianProcess{M,D,1,T,CartesianGrid{D,T},S}) where {M,D,T,S}
    for i in eachindex(gp.sim_prealloc.L,gp.sim_prealloc.Y)
        Z = complex(rand(Normal(0,1)), rand(Normal(0,1)))
        gp.sim_prealloc.Y[i] = gp.sim_prealloc.L[i] * Z
    end
    W = fft_array(gp.sim_prealloc.Y) / sqrt(length(gp.sim_prealloc.Y))
    return gp.mean[1] .+ real.(W[CartesianIndices(size(gp.mesh))])#, imag.(W[CartesianIndices(size(gp.mesh))])
end