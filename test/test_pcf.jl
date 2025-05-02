using ConservationLawsParticles
const PCF = PiecewiseConstantFunction

@testset "map" begin
    pcf = PCF([0., 1.], [0., 1., 2.])
    @test map(x -> x+1, pcf) == PCF([0., 1.], [1., 2., 3.])
end

@testset "abs" begin
    pcf = PCF([0., 1.], [-1., 0., 1.])
    @test abs(pcf) == PCF([0., 1.], [1., 0., 1.])
end

@testset "ops" begin
    f = PCF([   0.,    1., 2., 3.], [0., 1., 2., -1., 4.])
    g = PCF([-1., 0.5, 1.,  2.5],   [0., 1., -1., 0., 1.])
    @test f + g == PCF([-1., 0., 0.5, 1., 2., 2.5, 3.], [0., 1., 2., 0., 2., -1., 0., 5.])
    @test f - g == PCF([-1., 0., 0.5, 1., 2., 2.5, 3.], [0., -1., 0., 2., 2., -1., -2., 3.])
    @test f * g == PCF([-1., 0., 0.5, 1., 2., 2.5, 3.], [0., 0., 1., -1., 0., 0., -1., 4.])
    @test f ^ g == PCF([-1., 0., 0.5, 1., 2., 2.5, 3.], [1., 0., 1., 1., 1., 1., -1., 4.])
    # cannot use == directly for division because (f/g).values contains NaN
    f_g = f / g 
    @test f_g.cuts == [-1., 0., 0.5, 1., 2., 2.5, 3.]
    @test isnan(f_g.values[1])
    @test f_g.values[2:end] == [0., 1., -1., Inf, -Inf, -1., 4.]
end

@testset "lpnorm" begin
    pcf = PCF([0., 1., 11.], [0., 2., 3., 0.])
    @test lpnorm(pcf, 1) == 2 + 30
    @test lpnorm(pcf, 2) == sqrt(2^2 + 10*3^2)
    @test lpnorm(pcf, Inf) == 3

    pcf = PCF([0., 1.], [1., 1., 1.])
    # @test lpnorm(pcf, 1) == Inf
    @test lpnorm(pcf, Inf) == 1
end
