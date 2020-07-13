using Particles
using RecursiveArrayTools
using Test

@testset "densities" begin
    include("test_densities.jl")
end

@testset "velocities" begin
    include("test_velocities.jl")
end

@testset "make_velocities" begin
    V(x) = - x^3 + 0.2sin(12x)
    Wprime(r) = - 5 * sign(r) / (abs(r) + 1)
    mobility(rho) = max(1 - rho, 0)
    vel = make_velocity(V, Wprime, mobility)
    vels = make_velocities((V,), ((Wprime,),), (mobility,))
    x = ArrayPartition(sort(randn(500)))
    dx_vel = zero(x)
    dx_vels = zero(x)
    vel(dx_vel, x, nothing, 0.0)
    vels(dx_vels, x, nothing, 0.0)
    @test dx_vel == dx_vels
end

@testset "utils" begin
    @testset "empty_like" begin
        @test length(@test_deprecated empty_like([1,2,3])) == 3
        @test typeof(@test_deprecated empty_like([1,2,3])) == Vector{Int}
    end
end
