using ConservationLawsParticles
using RecursiveArrayTools, DifferentialEquations, Plots

# model
V(t, x) = 0.0 #2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 #x^2/2 - abs(x)
W′(t, x) = 0.0 #x - sign(x)
mob(ρ) = 0.0 # max((ρ-1)^2 * (2ρ+1), 0)
smodel = DiffusiveSampledModel((V,), ((W′,),), (mob,), (Diffusion(1),))
# imodel = DiffusiveIntegratedModel((V,), ((W,),), (mob,))

# initial condition
x0 = ArrayPartition(gaussian_particles(3, 201))

# time span
tspan = (0., 15.)

# ODE system for the particles
dprob = ODEProblem(velocities_diff!, x0, tspan, smodel)

abstol = reltol = 2e-8

@time dsol = solve(dprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(dsol, vars=1:5:201; color=:green)
savefig("diffusion1.png")

# plot an animation of the density
anim = @animate for t in range(0.0, 7.0; step=1/24)
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-10,10), ylims=(0,1/sqrt(pi)))
    densityplot!(dsol(t); color=:green, label="particles")
    plot!(x -> exp(-x^2/(1+4t))/sqrt(pi*(1+4t)); color=:orange, label="analytic")
end
gif(anim, "diffusion1.gif")
