export AbstractModel, SampledModel, IntegratedModel, DiffusiveSampledModel, DiffusiveIntegratedModel, HyperbolicModel, ParabolicModel
export num_species, external_velocity, interaction, mobility, diffusion, eachspecies, species
export SM, IM, DSM, DIM, HM, PM


"""
$(TYPEDEF)

Abstract type representing a particle model.
"""
abstract type AbstractModel{S} end


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
    SampledModel((V₁, ...), ((W′₁₁, ...), ...), (mob₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- sampled interactions `W′ᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`IntegratedModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = SampledModel(
           (V, V),
           ((Wprime_attr, Wprime_rep),
            (Wprime_rep, Wprime_attr)),
           (mobρ, mobσ));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities(x, model, 0.) .|> v -> round(v; digits=2)
([6.29, 0.25, -0.72, -7.23], [22.41, 1.12, 1.52, -7.87, -54.54])
```
"""
struct SampledModel{
    S,
    TVs         <: Tuple{Vararg{Any,S}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,S}},S}},
    Tmobilities <: Tuple{Vararg{Any,S}},
} <: AbstractModel{S}
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
end

SampledModel(V, Wprime, mob) = SampledModel((V,), ((Wprime,),), (mob,))


"""
    IntegratedModel((V₁, ...), ((W₁₁, ...), ...), (mob₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- integrated interactions `Wᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`SampledModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = IntegratedModel(
           (V, V),
           ((W_attr, W_rep),
            (W_rep, W_attr)),
           (mobρ, mobσ));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities(x, model, 0.) .|> v -> round(v; digits=2)
([6.62, 0.31, -0.43, -6.91], [23.26, 0.9, 0.71, -8.77, -55.23])
```
"""
struct IntegratedModel{
    S,
    TVs         <: Tuple{Vararg{Any,S}},
    TWs         <: Tuple{Vararg{Tuple{Vararg{Any,S}},S}},
    Tmobilities <: Tuple{Vararg{Any,S}},
} <: AbstractModel{S}
    Vs::TVs
    Ws::TWs
    mobilities::Tmobilities
end

IntegratedModel(V, W, mob) = IntegratedModel((V,), ((W,),), (mob,))


num_species(mod::SampledModel{S}) where S = S
num_species(mod::IntegratedModel{S}) where S = S

num_species(mod::Type{<:SampledModel{S}}) where S = S
num_species(mod::Type{<:IntegratedModel{S}}) where S = S

external_velocity(mod::SampledModel, i::Integer) = mod.Vs[i]
external_velocity(mod::IntegratedModel, i::Integer) = mod.Vs[i]

interaction(mod::SampledModel, i::Integer, j::Integer) = mod.Wprimes[i][j]
interaction(mod::IntegratedModel, i::Integer, j::Integer) = mod.Ws[i][j]

mobility(mod::SampledModel, i::Integer) = mod.mobilities[i]
mobility(mod::IntegratedModel, i::Integer) = mod.mobilities[i]

diffusion(mod::SampledModel, i::Integer) = nothing
diffusion(mod::IntegratedModel, i::Integer) = nothing


"""
    DiffusiveSampledModel((V₁, ...), ((W′₁₁, ...), ...), (mob₁, ...), (D₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- sampled interactions `W′ᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`,
- diffusions `Dᵢ`.

See also [`DiffusiveIntegratedModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = DiffusiveSampledModel(
           (V, V),
           ((Wprime_attr, Wprime_rep),
            (Wprime_rep, Wprime_attr)),
           (mobρ, mobσ),
           (Diffusion(1), nothing));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities_diff(x, model, 0.) .|> v -> round(v; digits=2)
([5.7, -0.81, 0.34, -6.64], [22.41, 1.12, 1.52, -7.87, -54.54])
```
"""
struct DiffusiveSampledModel{
    S,
    TVs         <: Tuple{Vararg{Any,S}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,S}},S}},
    Tmobilities <: Tuple{Vararg{Any,S}},
    Tdiffusions <: Tuple{Vararg{Any,S}}
} <: AbstractModel{S}
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
    diffusions::Tdiffusions
end

"""
    DiffusiveIntegratedModel((V₁, ...), ((W₁₁, ...), ...), (mob₁, ...), (D₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- integrated interactions `Wᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`,
- diffusions `Dᵢ`.

See also [`DiffusiveSampledModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = DiffusiveIntegratedModel(
           (V, V),
           ((W_attr, W_rep),
            (W_rep, W_attr)),
           (mobρ, mobσ),
           (Diffusion(1), nothing));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> velocities_diff(x, model, 0.) .|> v -> round(v; digits=2)
([6.03, -0.76, 0.63, -6.33], [23.26, 0.9, 0.71, -8.77, -55.23])
```
"""
struct DiffusiveIntegratedModel{
    S,
    TVs         <: Tuple{Vararg{Any,S}},
    TWs         <: Tuple{Vararg{Tuple{Vararg{Any,S}},S}},
    Tmobilities <: Tuple{Vararg{Any,S}},
    Tdiffusions <: Tuple{Vararg{Any,S}}
} <: AbstractModel{S}
    Vs::TVs
    Ws::TWs
    mobilities::Tmobilities
    diffusions::Tdiffusions
end


num_species(mod::DiffusiveSampledModel{S}) where S = S
num_species(mod::DiffusiveIntegratedModel{S}) where S = S

num_species(mod::Type{<:DiffusiveSampledModel{S}}) where S = S
num_species(mod::Type{<:DiffusiveIntegratedModel{S}}) where S = S

external_velocity(mod::DiffusiveSampledModel, i::Integer) = mod.Vs[i]
external_velocity(mod::DiffusiveIntegratedModel, i::Integer) = mod.Vs[i]

interaction(mod::DiffusiveSampledModel, i::Integer, j::Integer) = mod.Wprimes[i][j]
interaction(mod::DiffusiveIntegratedModel, i::Integer, j::Integer) = mod.Ws[i][j]

mobility(mod::DiffusiveSampledModel, i::Integer) = mod.mobilities[i]
mobility(mod::DiffusiveIntegratedModel, i::Integer) = mod.mobilities[i]

diffusion(mod::DiffusiveSampledModel, i::Integer) = mod.diffusions[i]
diffusion(mod::DiffusiveIntegratedModel, i::Integer) = mod.diffusions[i]


"""
    HyperbolicModel((V₁, ...), ((I₁₁, ...), ...), (mob₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- integrated interactions `Iᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`.

See also [`ParabolicModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = HyperbolicModel(
           (V, V),
           ((SampledInteraction(Wprime_attr), IntegratedInteraction(W_rep)),
            (IntegratedInteraction(W_rep), SampledInteraction(Wprime_attr))),
           (mobρ, mobσ));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> abstract_velocities(x, model, 0.) .|> v -> round(v; digits=2)
([6.27, 0.26, -0.38, -6.56], [22.92, 0.82, 0.71, -8.68, -54.89])
```
"""
struct HyperbolicModel{
    S,
    TVs <: NTuple{S,Any},
    TIs <: NTuple{S,NTuple{S,Any}},
    TMs <: NTuple{S,Any}
} <: AbstractModel{S}
    vels::TVs
    ints::TIs
    mobs::TMs
end


"""
    ParabolicModel((V₁, ...), ((I₁₁, ...), ...), (mob₁, ...), (D₁, ...))

Represents a particles system with:

- external velocities `Vᵢ`,
- interactions `Iᵢⱼ` (this is the effect of the species `j` on the species `i`),
- mobilities `mobᵢ`,
- diffusions `Dᵢ`.

See also [`HyperbolicModel`](@ref).

# Examples

```jldoctest
julia> using ConservationLawsParticles.Examples, RecursiveArrayTools

julia> model = ParabolicModel(
           (V, V),
           ((SampledInteraction(Wprime_attr), IntegratedInteraction(W_rep)),
            (IntegratedInteraction(W_rep), SampledInteraction(Wprime_attr))),
           (mobρ, mobσ),
           (Diffusion(1), nothing));

julia> x = ArrayPartition(gaussian_particles(2, 4), collect(range(-3, 4, length=5)));

julia> round.(x; digits=2)
([-2.0, -0.3, 0.3, 2.0], [-3.0, -1.25, 0.5, 2.25, 4.0])

julia> abstract_velocities(x, model, 0.) .|> v -> round(v; digits=2)
([5.68, -0.8, 0.68, -5.97], [22.92, 0.82, 0.71, -8.68, -54.89])
```
"""
struct ParabolicModel{
    S,
    TVs <: NTuple{S,Any},
    TIs <: NTuple{S,NTuple{S,Any}},
    TMs <: NTuple{S,Any},
    TDs <: NTuple{S,Any}
} <: AbstractModel{S}
    vels::TVs
    ints::TIs
    mobs::TMs
    difs::TDs
end


num_species(mod::HyperbolicModel{S}) where S = S
num_species(mod::ParabolicModel{S}) where S = S

num_species(mod::Type{<:HyperbolicModel{S}}) where S = S
num_species(mod::Type{<:ParabolicModel{S}}) where S = S

external_velocity(mod::HyperbolicModel, i::Integer) = mod.vels[i]
external_velocity(mod::ParabolicModel, i::Integer) = mod.vels[i]

interaction(mod::HyperbolicModel, i::Integer, j::Integer) = mod.ints[i][j]
interaction(mod::ParabolicModel, i::Integer, j::Integer) = mod.ints[i][j]

mobility(mod::HyperbolicModel, i::Integer) = mod.mobs[i]
mobility(mod::ParabolicModel, i::Integer) = mod.mobs[i]

diffusion(mod::HyperbolicModel, i::Integer) = nothing
diffusion(mod::ParabolicModel, i::Integer) = mod.difs[i]


const SM = SampledModel
const IM = IntegratedModel
const DSM = DiffusiveSampledModel
const DIM = DiffusiveIntegratedModel
const HM = HyperbolicModel
const PM = ParabolicModel
