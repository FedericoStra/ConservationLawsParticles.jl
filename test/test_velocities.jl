using .ConservationLawsParticles.Examples
using .ConservationLawsParticles: make_velocity, make_velocities

@testset "sampled" begin
@testset "1-S" begin
    vel = make_velocity(V, Wprime_attr, mob)
    vels = make_velocities((V,), ((Wprime_attr,),), (mob,))
    model = SampledModel((V,), ((Wprime_attr,),), (mob,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_vel = zero(x)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        vel(dx_vel, x, nothing, 0.0)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_vel == dx_vels
        @test dx_vel == dx_par
        @test dx_vel == dx_gen
        @test dx_vel == dx_abs_par
        @test dx_vel == dx_abs_gen
    end
end

@testset "2-S" begin
    model = SampledModel((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    vels = make_velocities((V, V2), ((Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_vels == dx_par
        @test dx_vels == dx_gen
        @test dx_vels == dx_abs_par
        @test dx_vels == dx_abs_gen
    end
end

@testset "3-S" begin
    mob(a, b, c) = a + b*c
    model = SampledModel((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mob, mob, mob))
    vels = make_velocities((V, V, V2), ((Wprime_rep, Wprime_attr, Wprime_attr), (Wprime_attr, Wprime_rep, Wprime_attr), (Wprime_attr, Wprime_attr, Wprime_rep)), (mob, mob, mob))
    @testset for lengths in [(2,2,2), (2,3,5), (3,5,7), (7,10,15), (50,75,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_vels = zero(x)
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        vels(dx_vels, x, nothing, 0.0)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_vels == dx_par
        @test dx_vels == dx_gen
        @test dx_vels == dx_abs_par
        @test dx_vels == dx_abs_gen
    end
end
end # sampled

@testset "integrated" begin
@testset "1-S" begin
    model = IntegratedModel((V,), ((W_attr,),), (mob,))
    @testset for len in [2, 5, 10, 50, 100]
        x = ArrayPartition(sort(randn(100)))
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_par == dx_gen
        @test dx_par == dx_abs_par
        @test dx_par == dx_abs_gen
    end
end

@testset "2-S" begin
    model = IntegratedModel((V, V2), ((W_rep, W_attr), (W_attr, W_rep)), (mobρ, mobσ))
    @testset for lengths in [(2,2), (2,3), (3,5), (7,10), (50,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_par == dx_gen
        @test dx_par == dx_abs_par
        @test dx_par == dx_abs_gen
    end
end

@testset "3-S" begin
    mob(a, b, c) = a + b*c
    model = IntegratedModel((V, V, V2), ((W_rep, W_attr, W_attr), (W_attr, W_rep, W_attr), (W_attr, W_attr, W_rep)), (mob, mob, mob))
    @testset for lengths in [(2,2,2), (2,3,5), (3,5,7), (7,10,15), (50,75,100)]
        x = ArrayPartition((sort(randn(len)) for len in lengths)...)
        dx_par = zero(x)
        dx_gen = zero(x)
        dx_abs_par = zero(x)
        dx_abs_gen = zero(x)
        velocities!(dx_par, x, model, 0.0)
        velocities_gen!(dx_gen, x, model, 0.0)
        abstract_velocities!(dx_abs_par, x, model, 0.0)
        abstract_velocities_gen!(dx_abs_gen, x, model, 0.0)
        @test dx_par == dx_gen
        @test dx_par == dx_abs_par
        @test dx_par == dx_abs_gen
    end
end
end # integrated

@testset "sampled-integrated" begin
@testset "1-S" begin
    smodel = SampledModel((V,), ((Wprime_attr,),), (mob,))
    imodel = IntegratedModel((V,), ((W_attr,),), (mob,))
    lengths = (200)
    x = ArrayPartition((sort(gaussian_particles(2, len)) for len in lengths)...)
    dx_s = zero(x)
    dx_i = zero(x)
    velocities_gen!(dx_s, x, smodel, 0.0)
    velocities_gen!(dx_i, x, imodel, 0.0)
    @test maximum(abs.(dx_s - dx_i)) < 0.01
end

@testset "2-S" begin
    smodel = SampledModel((V,V), ((Wprime_attr,Wprime_rep),(Wprime_rep,Wprime_attr)), (mobρ,mobσ))
    imodel = IntegratedModel((V,V), ((W_attr,W_rep),(W_rep,W_attr)), (mobρ,mobσ))
    lengths = (200, 200)
    x = ArrayPartition((sort(gaussian_particles(2, len)) for len in lengths)...)
    x.x[2] .+= 1
    dx_s = zero(x)
    dx_i = zero(x)
    velocities_gen!(dx_s, x, smodel, 0.0)
    velocities_gen!(dx_i, x, imodel, 0.0)
    @test maximum(abs.(dx_s - dx_i)) < 0.022
end
end # sampled-integrated
