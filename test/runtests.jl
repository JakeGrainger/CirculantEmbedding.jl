using CirculantEmbedding
import CirculantEmbedding.StaticArrays: SVector
import CirculantEmbedding.Meshes: CartesianGrid
using Test

@testset "CirculantEmbedding.jl" begin
    @testset "IndependentFields" begin
        grid = CartesianGrid((10,10), (0.,0.), (1.,1.))
        grid_alt = CartesianGrid((10,10), (0.,0.1), (1.,1.))
        @test_throws AssertionError IndependentFields((GaussianProcess(0.0, Matern(1., 1, 2, 2), grid, pad=2), GaussianProcess(0.0, Matern(1., 1, 2, 2), grid_alt, pad=2)))
        gp = IndependentFields(ntuple(d->GaussianProcess(0.0, Matern(1., 1, 2, 2), grid, pad=2), Val{2}()))
        @test rand(gp) isa Matrix{SVector{2,Float64}}
        @test getmesh(gp) == grid
    end
end
