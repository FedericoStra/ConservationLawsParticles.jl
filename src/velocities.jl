"""
    make_velocity(V::Function, Wprime::Function, mobility::Function)

Creates a function `velocity(dx, x, p, t)` that computes the velocity of the particles
under the influence of an external force `V`, a mutual interaction `Wprime` and the
congestion given by `mobility`.
"""
function make_velocity(V::Function, Wprime::Function, mobility::Function)
    function velocity(dx, x, p, t)
        R = pwc_density(x)
        len = length(x)
        for i in 1:len
            v::eltype(x) = 0
            for j in 1:i-1
                v += Wprime(x[j] - x[i])
            end
            for j in i+1:len
                v += Wprime(x[j] - x[i])
            end
            v /= len
            v += V(x[i])
            if v < 0
                mob = mobility(R[i])
            else
                mob = mobility(R[i+1])
            end
            dx[i] = v * mob
        end
    end
    velocity
end

export make_velocities_
function make_velocities_(
        Vs::Tuple{Vararg{Any,N}},
        Wprimes::Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
        mobilities::Tuple{Vararg{Any,N}}
        ) where N
    function velocities(
            dx::ArrayPartition{F, T},
            x::ArrayPartition{F, T},
            p,
            t
            ) where F where T<:Tuple{Vararg{AbstractVector{<:Real}}}
        dens = pwc_densities(x.x...)
        for spec in 1:N
            for i in 1:length(x.x[spec])
                v::eltype(dx.x[1]) = 0
                # self-interaction with `spec`
                local Wprime = Wprimes[spec][spec]
                for j in 1:i-1
                    v += Wprime(x.x[spec][j] - x.x[spec][i])
                end
                for j in i+1:length(x.x[spec])
                    v += Wprime(x.x[spec][j] - x.x[spec][i])
                end
                v /= length(x.x[spec])
                # interaction with `other < spec`
                for other in 1:spec-1
                    local Wprime = Wprimes[spec][other]
                    local w::eltype(dx.x[1]) = 0
                    for j in 1:length(x.x[other])
                        w += Wprime(x.x[other][j] - x.x[spec][i])
                    end
                    v += w / length(x.x[other])
                end
                # interaction with `other > spec`
                for other in spec+1:N
                    local Wprime = Wprimes[spec][other]
                    local w::eltype(dx.x[1]) = 0
                    for j in 1:length(x.x[other])
                        w += Wprime(x.x[other][j] - x.x[spec][i])
                    end
                    v += w / length(x.x[other])
                end
                # external velocity
                v += Vs[spec](x.x[spec][i])
                if v < 0
                    mob = mobilities[spec](dens[spec][:, 1, i]...)
                else
                    mob = mobilities[spec](dens[spec][:, 2, i]...)
                end
                dx.x[spec][i] = v * mob
            end
        end
    end
    velocities
end

export total_interaction, total_interaction_

function total_interaction(Wprime, ys::AbstractVector{<:Real}, x::Real)
    sum(Wprime(y - x) for y in ys) / length(ys)
end

function total_interaction_(x::Real; Wprime, particles::AbstractVector{<:Real})
    sum(Wprime(p - x) for p in particles) / length(particles)
end

function total_interaction_(Wprime, ys::AbstractVector{<:Real}, x::Real)
    T = promote_type(eltype(ys), typeof(x))
    w::T = 0
    for y in ys
        w += Wprime(y - x)
    end
    w / length(ys)
end

export integrated_total_interaction
function integrated_total_interaction(W, dens::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, x::Real)
    sum(i -> (dens[i] - dens[i+1]) * W(ys[i] - x), eachindex(ys))
end

function integrated_total_interaction(W, s::Integer, dens::AbstractArray{<:Real,3}, ys::AbstractVector{<:Real}, x::Real)
    sum(i -> (dens[s, 1, i] - dens[s, 2, i]) * W(ys[i] - x), eachindex(ys))
end

export fast_integrated_total_interaction
function fast_integrated_total_interaction(W, dens_diff::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, x::Real)
    sum(i -> dens_diff[i] * W(ys[i] - x), eachindex(ys))
end

function make_velocities(
        Vs::Tuple{Vararg{Any,N}},
        Wprimes::Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
        mobilities::Tuple{Vararg{Any,N}}
        ) where N
    function velocities(
            dx::ArrayPartition{F, T},
            x::ArrayPartition{F, T},
            p,
            t
            ) where F where T<:Tuple{Vararg{AbstractVector{<:Real}}}
        dens = pwc_densities(x.x...)
        for spec in 1:N
            for i in 1:length(x.x[spec])
                v::F = Vs[spec](x.x[spec][i])
                for other in 1:N
                    v += total_interaction(Wprimes[spec][other], x.x[other], x.x[spec][i])
                end
                if v < 0
                    mob = mobilities[spec](dens[spec][:, 1, i]...)
                else
                    mob = mobilities[spec](dens[spec][:, 2, i]...)
                end
                dx.x[spec][i] = v * mob
            end
        end
    end
    velocities
end

export Model
mutable struct Model{
    N,
    TVs         <: Tuple{Vararg{Any,N}},
    TWprimes    <: Tuple{Vararg{Tuple{Vararg{Any,N}},N}},
    Tmobilities <: Tuple{Vararg{Any,N}},
}
    Vs::TVs
    Wprimes::TWprimes
    mobilities::Tmobilities
end

export param_velocities
function param_velocities(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    dens = pwc_densities(x.x...)
    for spec in 1:N
        for i in 1:length(x.x[spec])
            v::F = p.Vs[spec](x.x[spec][i])
            for other in 1:N
                v += total_interaction(p.Wprimes[spec][other], x.x[other], x.x[spec][i])
            end
            if v < 0
                mob = p.mobilities[spec](dens[spec][:, 1, i]...)
            else
                mob = p.mobilities[spec](dens[spec][:, 2, i]...)
            end
            dx.x[spec][i] = v * mob
        end
    end
end

export param_velocities2
function param_velocities2(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    dens = pwc_densities(x.x...)
    for spec in 1:N
        dx.x[spec] .= p.Vs[spec].(x.x[spec])
        for other in 1:N
            for i in 1:length(x.x[spec])
                dx.x[spec][i] += total_interaction(p.Wprimes[spec][other], x.x[other], x.x[spec][i])
            end
        end
        for i in 1:length(x.x[spec])
            if dx.x[spec][i] < 0
                mob = p.mobilities[spec](dens[spec][:, 1, i]...)
            else
                mob = p.mobilities[spec](dens[spec][:, 2, i]...)
            end
            dx.x[spec][i] *= mob
        end
    end
end

export param_velocities3
function param_velocities3(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    dens = pwc_densities(x.x...)
    for spec in 1:N
        dx.x[spec] .= p.Vs[spec].(x.x[spec])
        for other in 1:N
        #     for i in 1:length(x.x[spec])
        #         dx.x[spec][i] += total_interaction_(x.x[spec][i]; Wprime=p.Wprimes[spec][other], particles=x.x[other])
        #     end
            dx.x[spec] .+= total_interaction_.(x.x[spec]; Wprime=p.Wprimes[spec][other], particles=x.x[other])
        end
        for i in 1:length(x.x[spec])
            if dx.x[spec][i] < 0
                mob = p.mobilities[spec](dens[spec][:, 1, i]...)
            else
                mob = p.mobilities[spec](dens[spec][:, 2, i]...)
            end
            dx.x[spec][i] *= mob
        end
    end
end

export gen_velocities
@generated function gen_velocities(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    init = :(dens = pwc_densities(x.x...))

    loop = Expr(:block, (quote
        for i in 1:length(x.x[$spec])
            v::F = p.Vs[$spec](x.x[$spec][i])

            $(Expr(:block, (quote
                v += total_interaction(p.Wprimes[$spec][$other], x.x[$other], x.x[$spec][i])
            end for other in 1:N)...))

            if v < 0
                mob = p.mobilities[$spec](dens[$spec][:, 1, i]...)
            else
                mob = p.mobilities[$spec](dens[$spec][:, 2, i]...)
            end
            dx.x[$spec][i] = v * mob
        end
    end for spec in 1:N)...)

    return :($init; $loop)
end

export gen_velocities2
@generated function gen_velocities2(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    init = :(dens = pwc_densities(x.x...))

    loop = Expr(:block, (quote
        dx.x[$spec] .= p.Vs[$spec].(x.x[$spec])
        $(Expr(:block, (quote
            for i in 1:length(x.x[$spec])
                dx.x[$spec][i] += total_interaction(p.Wprimes[$spec][$other], x.x[$other], x.x[$spec][i])
            end
        end for other in 1:N)...))
        for i in 1:length(x.x[$spec])
            if dx.x[$spec][i] < 0
                mob = p.mobilities[$spec](dens[$spec][:, 1, i]...)
            else
                mob = p.mobilities[$spec](dens[$spec][:, 2, i]...)
            end
            dx.x[$spec][i] *= mob
        end
    end for spec in 1:N)...)

    return :($init; $loop)
end

export gen_velocities3
@generated function gen_velocities3(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    init = :(dens = pwc_densities(x.x...))

    loop = Expr(:block, (quote
        dx.x[$spec] .= p.Vs[$spec].(x.x[$spec])
        $(Expr(:block, (quote
            for i in 1:length(x.x[$spec])
                dx.x[$spec][i] += total_interaction(p.Wprimes[$spec][$other], x.x[$other], x.x[$spec][i])
            end
        end for other in 1:N)...))
        for i in 1:length(x.x[$spec])
            d = dens[$spec]
            if dx.x[$spec][i] < 0
                mob = $(Expr(:call,
                    :(p.mobilities[$spec]),
                    (:(d[$j, 1, i]) for j in 1:N)...
                ))
            else
                mob = $(Expr(:call,
                    :(p.mobilities[$spec]),
                    (:(d[$j, 2, i]) for j in 1:N)...
                ))
            end
            dx.x[$spec][i] *= mob
        end
    end for spec in 1:N)...)

    return :($init; $loop)
end

export int_velocities
@generated function int_velocities(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::Model{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    init = quote
        dens = pwc_densities(x.x...)
        dens_diff = Vector{Vector{F}}(undef, N)
        for s in 1:N
            dens_diff[s] = dens[s][s, 1, :] - dens[s][s, 2, :]
        end
    end

    loop = Expr(:block, (quote
        dx.x[$spec] .= p.Vs[$spec].(x.x[$spec])
        $(Expr(:block, (quote
            for i in 1:length(x.x[$spec])
                dx.x[$spec][i] += fast_integrated_total_interaction(p.Wprimes[$spec][$other], dens_diff[$other], x.x[$other], x.x[$spec][i])
            end
        end for other in 1:N)...))
        for i in 1:length(x.x[$spec])
            d = dens[$spec]
            if dx.x[$spec][i] < 0
                mob = $(Expr(:call,
                    :(p.mobilities[$spec]),
                    (:(d[$j, 1, i]) for j in 1:N)...
                ))
            else
                mob = $(Expr(:call,
                    :(p.mobilities[$spec]),
                    (:(d[$j, 2, i]) for j in 1:N)...
                ))
            end
            dx.x[$spec][i] *= mob
        end
    end for spec in 1:N)...)

    return :($init; $loop)
end
