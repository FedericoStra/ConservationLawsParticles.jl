using ConservationLawsParticles
using RecursiveArrayTools, DifferentialEquations, Plots

# model
V(t, x) = t < 5 ? one(x) : -one(x)
Wprime(t, x) = sign(x) / (abs(x) + 1)
W(t, x) = log(abs(x) + 1)
mob(ρ) = max(1 - ρ, 0)
smodel = SampledModel((V,), ((Wprime,),), (mob,))
imodel = IntegratedModel((V,), ((W,),), (mob,))

# initial condition
x0 = ArrayPartition(vcat(range(0, 1, length=80)))

# time span
tspan = (0., 10.)

# ODE system for the particles
sprob = ODEProblem(velocities!, x0, tspan, smodel)
iprob = ODEProblem(velocities!, x0, tspan, imodel)

abstol = reltol = 1e-7;
@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);
@time isol = solve(iprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(title="Trajectories", xlabel="time", ylabel="position", legend=false)
plot!(ssol, vars=1:2:80; color=:blue)
plot!(isol, vars=1:2:80; color=:red)
savefig("trajectories.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; step=1/48)
    p = plot(title="Density", xlabel="position", ylabel="density",
        legend=false, xlims=(-2,4), ylims=(0,1))
    plot_density!(p, ssol(t); color=:blue)
    plot_density!(p, isol(t); color=:red)
    p
end
gif(anim, "density.gif")
