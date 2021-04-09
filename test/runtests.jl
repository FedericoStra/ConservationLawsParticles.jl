using ConservationLawsParticles
using RecursiveArrayTools
using Test

@testset "densities" begin
    include("test_densities.jl")
end

@testset "velocities" begin
    include("test_velocities.jl")
end

@testset "utils" begin
    @testset "empty_like" begin
        @test length(@test_deprecated empty_like([1,2,3])) == 3
        @test typeof(@test_deprecated empty_like([1,2,3])) == Vector{Int}
    end
end
