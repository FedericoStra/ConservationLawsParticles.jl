using ConservationLawsParticles
using RecursiveArrayTools, DifferentialEquations, Plots

# model
V(t, x) = 3 + (1+sin(3t)) * cos(4x)
Wprime(t, x) = 2 * sign(x) / (abs(x) + 1)
mob(ρ) = max(1 - ρ, 0)
model = SampledModel((V,), ((Wprime,),), (mob,))

# initial condition
x0 = ArrayPartition(gaussian_particles(2, 801))

# time span
tspan = (0., 5.)

# ODE system for the particles
prob = ODEProblem(velocities!, x0, tspan, model)

sol = solve(prob, BS5(); reltol=1e-8, abstol=1e-8)

# plot the particle trajectories
plot(sol, vars=1:20:801; legend=false, color=:blue,
    title="Trajectories", xlabel="time", ylabel="position")
savefig("traffic1.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; step=1/48)
    plot_density(sol(t); color=:blue, legend=false, xlims=(-2,12), ylims=(0,1),
        title="Density", xlabel="position", ylabel="density")
end
gif(anim, "traffic1.gif")
