@testset "1-S" begin
    vel = make_velocity(V, Wprime_attr, mobility)
    vels = make_velocities((V,), ((Wprime_attr,),), (mobility,))
    vels_ = make_velocities_((V,), ((Wprime_attr,),), (mobility,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_vel = zero(x)
        dx_vels = zero(x)
        dx_vels_ = zero(x)
        vel(dx_vel, x, nothing, 0.0)
        vels(dx_vels, x, nothing, 0.0)
        vels_(dx_vels_, x, nothing, 0.0)
        @test dx_vel == dx_vels == dx_vels_
    end
end

@testset "2-S" begin
    model = Model((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels = make_velocities((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels_ = make_velocities_((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        dx_vels_ = zero(x); vels_(dx_vels_, x, nothing, 0.0)
        @test_skip dx_vels == dx_vels_
        dx_p = zero(x); param_velocities(dx_p, x, model, 0.0)
        @test dx_vels == dx_p
        dx_p2 = zero(x); param_velocities2(dx_p2, x, model, 0.0)
        @test dx_vels == dx_p
        dx_g = zero(x); gen_velocities(dx_g, x, model, 0.0)
        @test dx_vels == dx_g
        dx_g2 = zero(x); gen_velocities2(dx_g2, x, model, 0.0)
        @test dx_vels == dx_g2
        dx_g3 = zero(x); gen_velocities3(dx_g3, x, model, 0.0)
        @test dx_vels == dx_g3
    end
end
