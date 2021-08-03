export sampled_interaction, integrated_interaction


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
