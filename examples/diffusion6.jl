using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots
using Random

# reference solution
boxsol(t, x; a, b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
refsol(t, x) = sum(boxsol(t, x; a=x0[i-1], b=x0[i]) for i in eachindex(x0)[2:end]) / (length(x0)-1)

# model
V(t, x) = 0.0 # 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = 0.0 # max((ρ-1)^2 * (2ρ+1), 0)
simple_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (SimpleDiffusion(1),))
mean_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MeanDiffusion(1),))
min_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MinDiffusion(1),))

# initial condition
Random.seed!(42)
n = 25
x0 = ArrayPartition(sort(vcat(randn(n) .- 2, .5randn(n) .+ 2)))

# time span
tspan = (0., 0.1)

# ODE system for the particles
simple_prob = ODEProblem(velocities_diff!, x0, tspan, simple_model)
mean_prob = ODEProblem(velocities_diff!, x0, tspan, mean_model)
min_prob = ODEProblem(velocities_diff!, x0, tspan, min_model)

abstol = reltol = 2e-7

@time simple_sol = solve(simple_prob, BS5(); abstol=abstol, reltol=reltol);
@time mean_sol = solve(mean_prob, BS5(); abstol=abstol, reltol=reltol);
@time min_sol = solve(min_prob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(title="Trajectories", xlabel="time", ylabel="position", legend=false)
plot!(simple_sol, vars=eachindex(x0)[1:1:end], color=1)
plot!(mean_sol, vars=eachindex(x0)[1:1:end], color=2)
plot!(min_sol, vars=eachindex(x0)[1:1:end], color=3)
savefig("diffusion6.png")

# plot an animation of the density
anim = @animate for t in range(sqrt.(tspan)...; length=100).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-5,3), ylims=(0,.8), legend=:topleft)
    plot!(x -> refsol(t, x); color=:black, label="analytic")
    densityplot!(simple_sol(t); color=1, label="simple")
    densityplot!(mean_sol(t); color=2, label="mean")
    densityplot!(min_sol(t); color=3, label="min")
end
gif(anim, "diffusion6.gif")

# plot an animation of the density
let tspan = (0., .003)
anim = @animate for t in range(sqrt.(tspan)...; length=100).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-5,3), ylims=(0,.8), legend=:topleft)
    plot!(x -> refsol(t, x); color=:black, label="analytic")
    densityplot!(simple_sol(t); color=1, label="simple")
    densityplot!(mean_sol(t); color=2, label="mean")
    densityplot!(min_sol(t); color=3, label="min")
end
gif(anim, "diffusion6.1.gif")
end
