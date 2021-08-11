using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots

using Random
Random.seed!(42)

# model
V(t, x) = -x^3
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = 1.0 # max((ρ-1)^2 * (2ρ+1), 0)
smodel = DiffusiveSampledModel((V,), ((W′,),), (mob,), (Diffusion(1),))
imodel = DiffusiveIntegratedModel((V,), ((W,),), (mob,), (Diffusion(1),))

# initial condition
n = 200
x0 = ArrayPartition(vcat(range(-1, 1; length=n+1)))
x0 = ArrayPartition(sort(randn(n+1)))

refsol(t, x) = sqrt(2)/gamma(1/4) * exp(-x^4/4)

# time span
tspan = (0., 2.0)

# ODE system for the particles
sprob = ODEProblem(velocities_diff!, x0, tspan, smodel)
# iprob = ODEProblem(velocities_diff!, x0, tspan, imodel)

abstol = reltol = 2e-8

@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(ssol, vars=eachindex(x0)[1:5:end]; color=:blue)
savefig("diffusion4.png")

# plot an animation of the density
anim = @animate for t in range(sqrt.(tspan)...; length=101).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-3,3), ylims=(0,1/sqrt(pi)))
    plot!(x -> refsol(t, x); color=:blue, label="asymptotic")
    densityplot!(ssol(t); color=:red, label="particles")
end
gif(anim, "diffusion4.gif")
