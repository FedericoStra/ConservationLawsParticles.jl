@testset "sampled" begin
@testset "1-S" begin
    vel = make_velocity(V, Wprime_attr, mobility)
    vels = make_velocities((V,), ((Wprime_attr,),), (mobility,))
    model = SampledInteraction((V,), ((Wprime_attr,),), (mobility,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_vel = zero(x)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        vel(dx_vel, x, nothing, 0.0)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_vel == dx_vels
        @test dx_vel == dx_par
        @test dx_vel == dx_gen
    end
end

@testset "2-S" begin
    model = SampledInteraction((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels = make_velocities((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_vels == dx_par
        @test dx_vels == dx_gen
    end
end

@testset "3-S" begin
    mobility(a, b, c) = a + b*c
    model = SampledInteraction((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mobility, mobility, mobility))
    vels = make_velocities((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mobility, mobility, mobility))
    @testset for lengths in [(2,2,2), (2,3,5), (3,5,7), (7,10,15), (50,75,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_vels == dx_par
        @test dx_vels == dx_gen
    end
end
end # sampled

@testset "integrated" begin
@testset "1-S" begin
    model = IntegratedInteraction((V,), ((W_attr,),), (mobility,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_par = zero(x)
        dx_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_par == dx_gen
    end
end

@testset "2-S" begin
    model = IntegratedInteraction((V, V2), ((W_rep, W_attr), (W_attr, W_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_par = zero(x)
        dx_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_par == dx_gen
    end
end

@testset "3-S" begin
    mobility(a, b, c) = a + b*c
    model = IntegratedInteraction((V, V, V2), ((W_rep, W_attr, W_attr), (W_attr, W_rep, W_attr), (W_attr, W_attr, W_rep)), (mobility, mobility, mobility))
    @testset for lengths in [(2,2,2), (2,3,5), (3,5,7), (7,10,15), (50,75,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_par = zero(x)
        dx_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        @test dx_par == dx_gen
    end
end
end # integrated
