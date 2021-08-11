using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots
using Random
Random.seed!(42)

# model
V(t, x) = 0.0 # 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = 0.0 # max((ρ-1)^2 * (2ρ+1), 0)
smodel = SampledModel((V,), ((W′,),), (mob,))
imodel = IntegratedModel((V,), ((W,),), (mob,))

# initial condition
x0 = ArrayPartition(sort(randn(50)))

boxsol(t, x; a, b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
refsol(t, x) = sum(boxsol(t, x; a=x0[i-1], b=x0[i])
    for i in eachindex(x0)[2:end]) / (length(x0)-1)

# time span
tspan = (0., 0.1)

# ODE system for the particles
sprob = ODEProblem(velocities_diff!, x0, tspan, smodel)
# iprob = ODEProblem(velocities_diff!, x0, tspan, imodel)

abstol = reltol = 2e-8

@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);
# @time isol = solve(iprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(ssol, vars=eachindex(x0)[1:5:end]; color=:blue)
savefig("diffusion5.png")

# plot an animation of the density
anim = @animate for t in range(sqrt.(tspan)...; length=101).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-5,3), ylims=(0,.8))
    plot!(x -> refsol(t, x); color=:blue, label="analytic")
    densityplot!(ssol(t); color=:red, label="particles")
end
gif(anim, "diffusion5.gif")

# plot an animation of the density
tspan = (0., .001)
anim = @animate for t in range(sqrt.(tspan)...; length=101).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-4,2), ylims=(0,1), legend=:topleft)
    plot!(x -> refsol(t, x); color=:blue, label="analytic")
    densityplot!(ssol(t); color=:red, label="particles")
end
gif(anim, "diffusion5.1.gif")
