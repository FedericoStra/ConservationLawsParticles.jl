using Particles
using Test

@testset "empty_like" begin
    @test length(@test_deprecated empty_like([1,2,3])) == 3
    @test typeof(@test_deprecated empty_like([1,2,3])) == Vector{Int}
end

@testset "pwc_density" begin
    @test pwc_density([0, 1, 3]) == [0.0, 0.5, 0.25, 0.0]
    @test pwc_density([0., 1.]) == [0.0, 1.0, 0.0]
    @test pwc_density([0.]) == [0.0, 0.0]
end
