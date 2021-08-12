using SpecialFunctions, RecursiveArrayTools, DifferentialEquations, Plots
using ConservationLawsParticles

# model
V1(t, x) = 1 + sin(x)/2
V2(t, x) = -1 - cos(x)/2
Wₐ′(t,x) = sign(x) / (abs(x) + 1) + x^3/20
Wᵣ(t, x) = 1 / (abs(x) + 1)
mob1(ρ, σ) = max(1 - ρ - σ/2, 0)
mob2(ρ, σ) = max(1 - ρ/2 - σ, 0)
attr = SampledInteraction(Wₐ′)
rep = IntegratedInteraction(Wᵣ)

model = ParabolicModel(
    (V1, V2),
    ((attr, rep), (rep, attr)),
    (mob1, mob2),
    (SimpleDiffusion(1/16), SimpleDiffusion(1/8)))

mmodel = ParabolicModel(
    (V1, V2),
    ((attr, rep), (rep, attr)),
    (mob1, mob2),
    (MinDiffusion(1/16), MinDiffusion(1/8)))

# initial condition
n = 80
x0 = ArrayPartition(
    vcat(range(-5, -3, length=n÷2), range(-2, -1, length=n÷2)),
    vcat(range(1, 2, length=n÷2), range(2.5, 3, length=n÷4), range(3.5, 4, length=n÷4)))

# time span
tspan = (0., 12.0)

# ODE system for the particles
prob = ODEProblem(abstract_velocities_gen!, x0, tspan, model)
mprob = ODEProblem(abstract_velocities_gen!, x0, tspan, mmodel)

# solve it
abstol = reltol = 1e-7
@time sol = solve(prob, BS5(); abstol=abstol, reltol=reltol);
@time msol = solve(mprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(sol, vars=1:1:n; color=1)
plot!(sol, vars=n+1:1:2n; color=2)
plot!(msol, vars=1:1:n; color=1, ls=:dot)
plot!(msol, vars=n+1:1:2n; color=2, ls=:dot)
savefig("multi-spec-2.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; length=100)
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-6,5), ylims=(0,1), legend=false)
    densityplot!(sol(t).x[1]; color=1)
    densityplot!(sol(t).x[2]; color=2)
    densityplot!(msol(t).x[1]; color=1, ls=:dot)
    densityplot!(msol(t).x[2]; color=2, ls=:dot)
end
gif(anim, "multi-spec-2.gif")
