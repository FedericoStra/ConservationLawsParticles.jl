export SampledInteraction, IntegratedInteraction
export sampled_interaction, integrated_interaction


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


function integrated_interaction(x::Real, W, ys::AbstractVector{<:Real}, dens_diff::AbstractVector{<:Real}=-diff(pwc_density(ys)))
    sum(i -> dens_diff[i] * W(ys[i] - x), eachindex(ys))
end

function integrated_interaction(t::Real, x::Real, W, ys::AbstractVector{<:Real}, dens_diff::AbstractVector{<:Real}=-diff(pwc_density(ys)))
    sum(i -> dens_diff[i] * W(t, ys[i] - x), eachindex(ys))
end
