export AbstractModel, SampledModel, IntegratedModel, DiffusiveSampledModel, DiffusiveIntegratedModel, HyperbolicModel
export num_species, external_velocity, interaction, mobility, diffusion, eachspecies, species


"""
$(TYPEDEF)

Abstract type representing a particle model.
"""
abstract type AbstractModel end


"""
$(TYPEDSIGNATURES)

Returns the number of species of the model.
"""
num_species(mod::AbstractModel) = error("unimplemented")
num_species(mod::Type{<:AbstractModel}) = error("unimplemented")

"""
$(TYPEDSIGNATURES)

Returns the external velocity field associated to the species `i`.
"""
external_velocity(mod::AbstractModel, i::Integer) = error("unimplemented")

"""
$(TYPEDSIGNATURES)

Returns the interaction exerted on the species `i` by the species `j`.
"""
interaction(mod::AbstractModel, i::Integer, j::Integer) = error("unimplemented")

"""
$(TYPEDSIGNATURES)

Returns the mobility associated to the species `i`.
"""
mobility(mod::AbstractModel, i::Integer) = error("unimplemented")

"""
$(TYPEDSIGNATURES)

Returns the diffusion associated to the species `i`.
"""
diffusion(mod::AbstractModel, i::Integer) = error("unimplemented")

"""
$(TYPEDSIGNATURES)

Returns an iterator over the indices of the species of the model.
"""
eachspecies(mod::AbstractModel) = 1:num_species(mod)
eachspecies(mod::Type{<:AbstractModel}) = 1:num_species(mod)

"""
$(TYPEDSIGNATURES)

Returns the particles associated to the species `i`.
"""
species(state, i::Integer) = error("unimplemented")
species(a::ArrayPartition, i::Integer) = a.x[i]

"""
    SampledModel((V₁, ...), ((W′₁₁, ...), ...), (mob₁, ...)

Represents a particles system with:

- external velocities `Vᵢ`,
- sampled interactions `W′ᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`IntegratedModel`](@ref).

# Examples

```jldoctest; setup = :(using ConservationLawsParticles)
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = SampledModel(
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
struct SampledModel{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
} <: AbstractModel
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
end

SampledModel(V, Wprime, mob) = SampledModel((V,), ((Wprime,),), (mob,))


"""
    IntegratedModel((V₁, ...), ((W₁₁, ...), ...), (mob₁, ...)

Represents a particles system with:

- external velocities `Vᵢ`,
- integrated interactions `Wᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`SampledModel`](@ref).

# Examples

```jldoctest; setup = :(using ConservationLawsParticles)
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = IntegratedModel(
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
struct IntegratedModel{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWs         <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
} <: AbstractModel
    Vs::TVs
    Ws::TWs
    mobilities::Tmobilities
end

IntegratedModel(V, W, mob) = IntegratedModel((V,), ((W,),), (mob,))


num_species(mod::SampledModel{N}) where N = N
num_species(mod::IntegratedModel{N}) where N = N

num_species(mod::Type{<:SampledModel{N}}) where N = N
num_species(mod::Type{<:IntegratedModel{N}}) where N = N

external_velocity(mod::SampledModel, i::Integer) = mod.Vs[i]
external_velocity(mod::IntegratedModel, i::Integer) = mod.Vs[i]

interaction(mod::SampledModel, i::Integer, j::Integer) = mod.Wprimes[i][j]
interaction(mod::IntegratedModel, i::Integer, j::Integer) = mod.Ws[i][j]

mobility(mod::SampledModel, i::Integer) = mod.mobilities[i]
mobility(mod::IntegratedModel, i::Integer) = mod.mobilities[i]

diffusion(mod::SampledModel, i::Integer) = nothing
diffusion(mod::IntegratedModel, i::Integer) = nothing


struct DiffusiveSampledModel{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
    Tdiffusions <: Tuple{Vararg{Any,N}}
} <: AbstractModel
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
    diffusions::Tdiffusions
end

struct DiffusiveIntegratedModel{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWs         <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
    Tdiffusions <: Tuple{Vararg{Any,N}}
} <: AbstractModel
    Vs::TVs
    Ws::TWs
    mobilities::Tmobilities
    diffusions::Tdiffusions
end


num_species(mod::DiffusiveSampledModel{N}) where N = N
num_species(mod::DiffusiveIntegratedModel{N}) where N = N

num_species(mod::Type{<:DiffusiveSampledModel{N}}) where N = N
num_species(mod::Type{<:DiffusiveIntegratedModel{N}}) where N = N

external_velocity(mod::DiffusiveSampledModel, i::Integer) = mod.Vs[i]
external_velocity(mod::DiffusiveIntegratedModel, i::Integer) = mod.Vs[i]

interaction(mod::DiffusiveSampledModel, i::Integer, j::Integer) = mod.Wprimes[i][j]
interaction(mod::DiffusiveIntegratedModel, i::Integer, j::Integer) = mod.Ws[i][j]

mobility(mod::DiffusiveSampledModel, i::Integer) = mod.mobilities[i]
mobility(mod::DiffusiveIntegratedModel, i::Integer) = mod.mobilities[i]

diffusion(mod::DiffusiveSampledModel, i::Integer) = mod.diffusions[i]
diffusion(mod::DiffusiveIntegratedModel, i::Integer) = mod.diffusions[i]


struct HyperbolicModel{
    N,
    TVs <: NTuple{N,Any},
    TIs <: NTuple{N,NTuple{N,Any}},
    TMs <: NTuple{N,Any}
} <: AbstractModel
    vels::TVs
    ints::TIs
    mobs::TMs
end


num_species(mod::HyperbolicModel{N}) where N = N

num_species(mod::Type{<:HyperbolicModel{N}}) where N = N

external_velocity(mod::HyperbolicModel, i::Integer) = mod.vels[i]

interaction(mod::HyperbolicModel, i::Integer, j::Integer) = mod.ints[i][j]

mobility(mod::HyperbolicModel, i::Integer) = mod.mobs[i]

diffusion(mod::HyperbolicModel, i::Integer) = nothing
