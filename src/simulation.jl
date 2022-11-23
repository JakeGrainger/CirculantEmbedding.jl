struct GaussianProcess{D,P,T<:Real,G<:Mesh{D,T},S}
    Γ::Kernel{D,P}
    mesh::G
    sim_prealloc::S
end
GaussianProcess(Γ::Kernel, mesh::Mesh; pad=0) = GaussianProcess(Γ, mesh, preallocate_gp_simulation(Γ, mesh, pad))

preallocate_gp_simulation(Γ::Kernel, mesh::Mesh, pad) = error("Simulation only currently implemented on a regular grid.")
preallocate_gp_simulation(Γ::Kernel{D,P}, mesh::CartesianGrid{D,T}, pad) where {D,P,T} = CirculantPrealloc(Γ, mesh, pad)
preallocate_gp_simulation(Γ::Kernel{D,P}, mesh::CartesianGrid{D2,T}, pad) where {D,D2,P,T} = error("Kernel has different dimension to grid.")

struct CirculantPrealloc{T,S,D}
    L::Array{T,D}
    Y::Array{S,D}
    function CirculantPrealloc(Γ::Kernel{D,P}, mesh::CartesianGrid{D,T}, pad=0) where {D,P,T}
        C = cov_at_ft(Γ, mesh, pad)
        A = fft_array(C)
        L = cholesky.(A)
        Y = Array{SVector{P,ComplexF64}}(undef, size(L))
        new{eltype(L), eltype(Y), D}(L,Y)
    end
    function CirculantPrealloc(Γ::Kernel{D,1}, mesh::CartesianGrid{D,T}, pad=0) where {D,T}
        C = cov_at_ft(Γ, mesh, pad)
        A = fft(C)
        L = sqrt.(A)
        Y = Array{ComplexF64}(undef, size(L))
        new{eltype(L), eltype(Y), D}(L,Y)
    end
end

function cov_at_ft(Γ::Kernel, mesh, pad=0)
    lags = compute_lags(mesh, pad)
    return SHermitianCompact.(cov.(Γ, lags)) # strictly not symmetric, but doesn't change the algorithm as result of fft will be hermitian symm.
end

const KernelSomeSdfOnly = Union{KernelSdfOnly, AdditiveKernel{<:KernelSdfOnly,<:Kernel}, AdditiveKernel{<:Kernel, <:KernelSdfOnly}, AdditiveKernel{<:KernelSdfOnly, <:KernelSdfOnly}}
function cov_at_ft(Γ::Union{KernelSdfOnly, AdditiveKernel}, mesh, pad)
    lags = compute_lags(mesh, pad)
    return SHermitianCompact.(approx_cov(Γ, lags))
end

function compute_lags(mesh::CartesianGrid{D,T}, pad) where {D,T}
    n, Δ = mesh.dims, mesh.spacing
    m = ntuple(d -> choose_circ_m(n[d], pad), Val{D}())
    Iterators.product(ntuple(d -> fftfreq(m[d],m[d]*Δ[d]), Val{D}())...)
end

function choose_circ_m(n::Int, pad)
    return 2n-2+pad
end

function rand(gp::GaussianProcess{D,P,T,CartesianGrid{D,T},S}) where {D,P,T,S}
    for i in eachindex(gp.sim_prealloc.L,gp.sim_prealloc.Y)
        Z = SVector{P,complex(T)}(complex(rand(Normal(0,1)), rand(Normal(0,1))) for i in 1:P)
        gp.sim_prealloc.Y[i] = gp.sim_prealloc.L[i].L * Z
    end
    W = fft_array(gp.sim_prealloc.Y) / sqrt(length(gp.sim_prealloc.Y))
    return real.(W[CartesianIndices(size(gp.mesh))])#, imag.(W[CartesianIndices(size(gp.mesh))])
end
function rand(gp::GaussianProcess{D,1,T,CartesianGrid{D,T},S}) where {D,T,S}
    for i in eachindex(gp.sim_prealloc.L,gp.sim_prealloc.Y)
        Z = complex(rand(Normal(0,1)), rand(Normal(0,1)))
        gp.sim_prealloc.Y[i] = gp.sim_prealloc.L[i] * Z
    end
    W = fft_array(gp.sim_prealloc.Y) / sqrt(length(gp.sim_prealloc.Y))
    return real.(W[CartesianIndices(size(gp.mesh))])#, imag.(W[CartesianIndices(size(gp.mesh))])
end