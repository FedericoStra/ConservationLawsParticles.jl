using Particles
using RecursiveArrayTools
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

@testset "pwc_densities" begin
    xs = ([0., 1., 2.], [0., 1., 2., 3.])
    @test pwc_density(xs...) == pwc_densities(xs...)
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
