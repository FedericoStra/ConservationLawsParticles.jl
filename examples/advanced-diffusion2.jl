using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots

# model
V(t, x) = 0.0 # 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = 0.0 # max((ρ-1)^2 * (2ρ+1), 0)
smodel = SampledModel((V,), ((W′,),), (mob,))
imodel = IntegratedModel((V,), ((W,),), (mob,))

# initial condition
n = 20
x0 = ArrayPartition(vcat(
    range(-3, -2; length=2n),
    range(-1, 0; length=n),
    range(1, 2; length=3n)
    ))

boxsol(t, x; a, b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
refsol(t, x) = (
        boxsol(t, x; a=-3, b=-2) * (2n-1) +
        boxsol(t, x; a=-1, b=0) * (n-1) +
        boxsol(t, x; a=1, b=2) * (3n-1) +
        boxsol(t, x; a=-2, b=-1) + boxsol(t, x; a=0, b=1)
    ) / (6n-1)

# time span
tspan = (0., 1.0)

# ODE system for the particles
sprob = ODEProblem(velocities_diff!, x0, tspan, smodel)
# iprob = ODEProblem(velocities_diff!, x0, tspan, imodel)

abstol = reltol = 2e-8

@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(ssol, vars=eachindex(x0)[1:5:end]; color=:blue)
savefig("advanced-diffusion2.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; length=101)
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-7,7), ylims=(0,1/sqrt(pi)))
    plot!(x -> refsol(t, x); color=:blue, label="analytic")
    densityplot!(ssol(t); color=:red, label="particles")
end
gif(anim, "advanced-diffusion2.gif")
