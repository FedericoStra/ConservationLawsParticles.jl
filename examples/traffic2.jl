using ConservationLawsParticles
using RecursiveArrayTools, DifferentialEquations, Plots

# model
V(t, x) = 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = x^2/2 - abs(x)
W′(t, x) = x - sign(x)
# W(t, x) = (a=abs(x); a*log(a+1/4))
# W′(t, x) = (a=abs(x); sign(x) * (log(a+1/4) + a/(a+1/4)))
mob(ρ) = max((ρ-1)^2 * (2ρ+1), 0)
smodel = SampledModel((V,), ((W′,),), (mob,))
imodel = IntegratedModel((V,), ((W,),), (mob,))

# initial condition
x0 = ArrayPartition(gaussian_particles(2, 401))

# time span
tspan = (0., 15.)

# ODE system for the particles
sprob = ODEProblem(velocities_gen!, x0, tspan, smodel)
iprob = ODEProblem(velocities_gen!, x0, tspan, imodel)

abstol = reltol = 2e-8

@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);
@time isol = solve(iprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(ssol, vars=1:20:401; color=:blue)
plot!(isol, vars=1:20:401; color=:red)
savefig("traffic2.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; step=1/24)
    p = plot(legend=false, xlims=(-2,20), ylims=(0,1),
        title="Density", xlabel="position", ylabel="density")
    x = ssol(t).x[1]
    min_dens = minimum(pwc_densities(x)[1][1,:,:], dims=1)[:];
    plot!(p, repeat(x, inner=(1,2))[1:20:401,:]',
        hcat(zero(min_dens), min_dens)[1:20:401,:]',
        color=:gray, opacity=.25)
    plot_density!(p, ssol(t); color=:blue)
    plot_density!(p, isol(t); color=:red)
end
gif(anim, "traffic2.gif")
