using ConservationLawsParticles

@testset "@time_independent" begin
    @time_independent V_time_indep(x, y) = x + y
    @test V_time_indep(2, 3) == 5
    @test V_time_indep(0, 2, 3) == 5
    @test V_time_indep(:time, 2, 3) == 5
    @test V_time_indep(nothing, 2, 3) == 5
end
