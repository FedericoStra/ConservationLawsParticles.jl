export SampledInteraction, IntegratedInteraction
export sampled_interaction, integrated_interaction

"""
    SampledInteraction((V₁, ...), ((W′₁₁, ...), ...), (mob₁, ...)

Represents a particles system with:

- external velocities `Vᵢ`,
- sampled interactions `W′ᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`IntegratedInteraction`](@ref).

# Examples

```jldoctest; setup = :(using ConservationLawsParticles)
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = SampledInteraction(
           (V, V),
           ((Wprime_attr, Wprime_rep),
            (Wprime_rep, Wprime_attr)),
           (mobρ, mobσ));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)))
([-2.000000000000003, -0.30305473018369145, 0.30305473018369145, 2.000000000000003], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities(x, model, 0.)
([6.291136247066298, 0.2466663161150116, -0.7218917091339228, -7.22630670503873], [22.405129478914613, 1.1249366684885518, 1.5188519354999799, -7.87111869358889, -54.536397957423915])
```
"""
mutable struct SampledInteraction{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
}
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
end


"""
    IntegratedInteraction((V₁, ...), ((W₁₁, ...), ...), (mob₁, ...)

Represents a particles system with:

- external velocities `Vᵢ`,
- integrated interactions `Wᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`SampledInteraction`](@ref).

# Examples

```jldoctest; setup = :(using ConservationLawsParticles)
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = IntegratedInteraction(
           (V, V),
           ((W_attr, W_rep),
            (W_rep, W_attr)),
           (mobρ, mobσ));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)))
([-2.000000000000003, -0.30305473018369145, 0.30305473018369145, 2.000000000000003], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities(x, model, 0.)
([6.621647425385332, 0.30545649966452776, -0.42671180044513507, -6.914302510772401], [23.261098611816987, 0.9023583690557855, 0.7055593913708742, -8.774526834146023, -55.23308442021724])
```
"""
mutable struct IntegratedInteraction{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWs         <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
}
    Vs::TVs
    Ws::TWs
    mobilities::Tmobilities
end


@doc raw"""
    sampled_interaction([t,] x, W′, ys)

Computes ``-(W' * \dot\rho)(t, x)``, which is the sampled approximation of ``-(W' * \rho)(t, x)``,
where ``\rho`` is the piecewise-constant density associated to the particles `ys`.

It `t` is omitted, then `W′(x)` is assumed independent of time.

See also [`integrated_interaction`](@ref).
"""
function sampled_interaction end

function sampled_interaction(x::Real, Wprime, ys::AbstractVector{<:Real})
    -sum(Wprime(x - y) for y in ys) / (length(ys) - 1)
end

function sampled_interaction(t::Real, x::Real, Wprime, ys::AbstractVector{<:Real})
    -sum(Wprime(t, x - y) for y in ys) / (length(ys) - 1)
end

function sampled_interaction(x::Real; Wprime, particles::AbstractVector{<:Real})
    -sum(Wprime(x - p) for p in particles) / (length(particles) - 1)
end

function sampled_interaction(t::Real, x::Real; Wprime, particles::AbstractVector{<:Real})
    -sum(Wprime(t, x - p) for p in particles) / (length(particles) - 1)
end

# function total_interaction__(Wprime, ys::AbstractVector{<:Real}, x::Real)
#     T = promote_type(eltype(ys), typeof(x))
#     w::T = 0
#     for y in ys
#         w += Wprime(y - x)
#     end
#     w / (length(ys) - 1)
# end


@doc raw"""
    integrated_interaction([t,] x, W, ys[, dens_diff])

Computes ``-(W' * \rho)(t, x) = -(W * \rho')(t,x)``, where ``\rho`` is the piecewise-constant
density associated to the particles `ys`.

It `t` is omitted, then `W(x)` is assumed independent of time.

!!! note

    To ensure the correctness of the computation, `dens_diff` must coincide with `diff(pwc_density(ys))`.
    It can be pre-computed and passed explicitly to allow reuse (as an optimization).

See also [`sampled_interaction`](@ref).
"""
function integrated_interaction end

function integrated_interaction(x::Real, W, ys::AbstractVector{<:Real}, dens_diff::AbstractVector{<:Real}=diff(pwc_density(ys)))
    -sum(i -> dens_diff[i] * W(x - ys[i]), eachindex(ys))
end

function integrated_interaction(t::Real, x::Real, W, ys::AbstractVector{<:Real}, dens_diff::AbstractVector{<:Real}=diff(pwc_density(ys)))
    -sum(i -> dens_diff[i] * W(t, x - ys[i]), eachindex(ys))
end
