@testset "1-S" begin
    vel = make_velocity(V, Wprime_attr, mobility)
    vels = make_velocities((V,), ((Wprime_attr,),), (mobility,))
    vels_ = make_velocities_((V,), ((Wprime_attr,),), (mobility,))
    model = Model((V,), ((Wprime_attr,),), (mobility,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_vel = zero(x)
        dx_vels = zero(x)
        dx_vels_ = zero(x)
        dx_param = zero(x)
        dx_param2 = zero(x)
        dx_gen = zero(x)
        dx_gen2 = zero(x)
        dx_gen3 = zero(x)
        vel(dx_vel, x, nothing, 0.0)
        vels(dx_vels, x, nothing, 0.0)
        vels_(dx_vels_, x, nothing, 0.0)
        param_velocities(dx_param, x, model, 0.0)
        param_velocities2(dx_param2, x, model, 0.0)
        gen_velocities(dx_gen, x, model, 0.0)
        gen_velocities2(dx_gen2, x, model, 0.0)
        gen_velocities3(dx_gen3, x, model, 0.0)
        @test dx_vel == dx_vels
        @test dx_vel == dx_vels_
        @test dx_vel == dx_param
        @test dx_vel == dx_param2
        @test dx_vel == dx_gen
        @test dx_vel == dx_gen2
        @test dx_vel == dx_gen3
    end
end

@testset "2-S" begin
    model = Model((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels = make_velocities((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels_ = make_velocities_((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_vels_ = zero(x)
        dx_param = zero(x)
        dx_param2 = zero(x)
        dx_gen = zero(x)
        dx_gen2 = zero(x)
        dx_gen3 = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        vels_(dx_vels_, x, nothing, 0.0)
        param_velocities(dx_param, x, model, 0.0)
        param_velocities2(dx_param2, x, model, 0.0)
        gen_velocities(dx_gen, x, model, 0.0)
        gen_velocities2(dx_gen2, x, model, 0.0)
        gen_velocities3(dx_gen3, x, model, 0.0)
        @test_skip dx_vels == dx_vels_
        @test dx_vels == dx_param
        @test dx_vels == dx_param2
        @test dx_vels == dx_gen
        @test dx_vels == dx_gen2
        @test dx_vels == dx_gen3
    end
end

@testset "3-S" begin
    mobility(a, b, c) = a + b*c
    model = Model((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mobility, mobility, mobility))
    vels = make_velocities((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mobility, mobility, mobility))
    vels_ = make_velocities_((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mobility, mobility, mobility))
    @testset for lengths in [(2,2,2), (2,3,5), (3,5,7), (7,10,15), (50,75,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_vels_ = zero(x)
        dx_param = zero(x)
        dx_param2 = zero(x)
        dx_gen = zero(x)
        dx_gen2 = zero(x)
        dx_gen3 = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        vels_(dx_vels_, x, nothing, 0.0)
        param_velocities(dx_param, x, model, 0.0)
        param_velocities2(dx_param2, x, model, 0.0)
        gen_velocities(dx_gen, x, model, 0.0)
        gen_velocities2(dx_gen2, x, model, 0.0)
        gen_velocities3(dx_gen3, x, model, 0.0)
        @test_skip dx_vels == dx_vels_
        @test dx_vels == dx_param
        @test dx_vels == dx_param2
        @test dx_vels == dx_gen
        @test dx_vels == dx_gen2
        @test dx_vels == dx_gen3
    end
end
