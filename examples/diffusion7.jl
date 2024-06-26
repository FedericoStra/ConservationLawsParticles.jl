using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots
using Random

# reference solution
boxsol(t, x; a, b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
# refsol(t, x) = sum(boxsol(t, x; a=x0[i-1], b=x0[i]) for i in eachindex(x0)[2:end]) / (length(x0)-1)
refsol(t, x) = sqrt(2)/gamma(1/4) * exp(-x^4/4)

# model
V(t, x) = -x^3 # 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = 1.0 # max((ρ-1)^2 * (2ρ+1), 0)
simple_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (SimpleDiffusion(1),))
mean_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MeanDiffusion(1),))
min_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MinDiffusion(1),))

# initial condition
Random.seed!(42)
n = 40
x0 = ArrayPartition(sort(vcat(randn(n) .- 1, .5randn(n) .+ 1)))

# time span
tspan = (0., 3.)

# ODE system for the particles
simple_prob = ODEProblem(velocities_diff!, x0, tspan, simple_model)
mean_prob = ODEProblem(velocities_diff!, x0, tspan, mean_model)
min_prob = ODEProblem(velocities_diff!, x0, tspan, min_model)

abstol = reltol = 2e-7

@time simple_sol = solve(simple_prob, BS5(); abstol=abstol, reltol=reltol);
@time mean_sol = solve(mean_prob, BS5(); abstol=.5abstol, reltol=.5reltol);
@time min_sol = solve(min_prob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(title="Trajectories", xlabel="time", ylabel="position", legend=false)
plot!(simple_sol, vars=eachindex(x0)[1:5:end], color=1)
plot!(mean_sol, vars=eachindex(x0)[1:5:end], color=2)
plot!(min_sol, vars=eachindex(x0)[1:5:end], color=3)
savefig("diffusion7.png")

# plot an animation of the density
anim = @animate for t in range(sqrt.(tspan)...; length=100).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-3,3), ylims=(0,.4), legend=:topleft, size=(1200, 800))
    plot!(x -> refsol(t, x); color=:black, label="analytic")
    densityplot!(simple_sol(t); color=1, label="simple")
    densityplot!(mean_sol(t); color=2, label="mean")
    densityplot!(min_sol(t); color=3, label="min")
end
gif(anim, "diffusion7.gif")

# plot an animation of the density
let tspan = (0., .1)
anim = @animate for t in range(sqrt.(tspan)...; length=100).^2
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-3,3), ylims=(0,.4), legend=:topleft)
    plot!(x -> refsol(t, x); color=:black, label="analytic")
    densityplot!(simple_sol(t); color=1, label="simple")
    densityplot!(mean_sol(t); color=2, label="mean")
    densityplot!(min_sol(t); color=3, label="min")
end
gif(anim, "diffusion7.1.gif")
end
