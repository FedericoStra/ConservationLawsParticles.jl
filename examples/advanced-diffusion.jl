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
n = 100
x0 = ArrayPartition(vcat(
    range(-2, -1; length=n+1),
    range(0, 2; length=n+1)
    ))

boxsol(t,x;a,b) = (erf((x-a)/2sqrt(t)) - erf((x-b)/2sqrt(t))) / (2*(b-a))
refsol(t, x) = boxsol(t,x; a=-2,b=-1) * n/(2n+1) +
    boxsol(t,x; a=0,b=2) * n/(2n+1) +
    boxsol(t,x; a=-1,b=0) * 1/(2n-1)

# time span
tspan = (0., 2.0)

# ODE system for the particles
sprob = ODEProblem(velocities_diff!, x0, tspan, smodel)
# iprob = ODEProblem(velocities_diff!, x0, tspan, imodel)

abstol = reltol = 2e-8

@time ssol = solve(sprob, BS5(); abstol=abstol, reltol=reltol);
# @time isol = solve(iprob, BS5(); abstol=abstol, reltol=reltol);

# plot the particle trajectories
plot(legend=false, title="Trajectories", xlabel="time", ylabel="position")
plot!(ssol, vars=1:5:(2n+2); color=:blue)
savefig("advanced-diffusion.png")

# plot an animation of the density
anim = @animate for t in range(tspan...; length=101)
    plot(title="Density", xlabel="position", ylabel="density",
        xlims=(-7,7), ylims=(0,1/sqrt(pi)))
    plot!(x -> refsol(t, x); color=:blue, label="analytic")
    densityplot!(ssol(t); color=:red, label="particles")
end
gif(anim, "advanced-diffusion.gif")
