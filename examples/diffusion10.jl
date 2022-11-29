using ConservationLawsParticles
using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots
using Random

# reference solution
boxsol(t, x; a, b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
# refsol(t, x) = sum(boxsol(t, x; a=x0[i-1], b=x0[i]) for i in eachindex(x0)[2:end]) / (length(x0)-1)
# refsol(t, x) = 1 / (pi * (1 + x^2))

# model
V(t, x) = sin(3t) - (x/1.75)^5 # -2x / (1 + x^2) # 2 + 0.5sin(2x) + 0.3sin(2sqrt(2)*x) + 0.2cos(2sqrt(7)*x)
W(t, x) = 0.0 # x^2/2 - abs(x)
W′(t, x) = 0.0 # x - sign(x)
mob(ρ) = max((ρ-1)^2 * (2ρ+1), 0)
function A(r)
    if r < 1/4
        r
    elseif r < 3/4
        1/4
    else
        r - 1/2
    end
end
simple_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (SimpleDiffusion(A),))
mean_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MeanDiffusion(A),))
min_model = DiffusiveSampledModel((V,), ((W′,),), (mob,), (MinDiffusion(A),))

# initial condition
Random.seed!(42)
n = 40
x0 = ArrayPartition(sort(vcat(.1randn(n) .- .5, .25randn(n) .+ .5)))

# time span
tspan = (0., 5.)

# ODE system for the particles
simple_prob = ODEProblem(velocities_diff!, x0, tspan, simple_model)
mean_prob = ODEProblem(velocities_diff!, x0, tspan, mean_model)
min_prob = ODEProblem(velocities_diff!, x0, tspan, min_model)

# solve it
abstol = reltol = 1e-7
@time simple_sol = solve(simple_prob, BS5(); abstol=abstol, reltol=reltol);
@time mean_sol = solve(mean_prob, BS5(); abstol=abstol, reltol=reltol);
@time min_sol = solve(min_prob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(title="Trajectories", xlabel="time", ylabel="position", legend=false, ylims=(-5,5))
plot!(simple_sol, vars=eachindex(x0)[1:2:end], color=1)
plot!(mean_sol, vars=eachindex(x0)[1:2:end], color=2)
plot!(min_sol, vars=eachindex(x0)[1:2:end], color=3)
savefig("diffusion10.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; length=100)
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-4,4), ylims=(0,1), legend=:topleft, size=(1200, 800))
    # plot!(x -> refsol(t, x); color=:black, label="analytic")
    densityplot!(simple_sol(t); color=1, label="simple")
    densityplot!(mean_sol(t); color=2, label="mean")
    densityplot!(min_sol(t); color=3, label="min")
end
gif(anim, "diffusion10.gif")
