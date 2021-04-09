using ConservationLawsParticles
using RecursiveArrayTools
using Test
using TestSetExtensions

@testset ExtendedTestSet "ConservationLawsParticles.jl" begin
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
        @testset "@time_independent" begin
            @time_independent V_time_indep(x, y) = x + y
            @test V_time_indep(2, 3) == 5
            @test V_time_indep(0, 2, 3) == 5
            @test V_time_indep(:time, 2, 3) == 5
            @test V_time_indep(nothing, 2, 3) == 5
        end
    end
end
