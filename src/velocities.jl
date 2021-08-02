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
            v /= len - 1
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
                v /= length(x.x[spec]) - 1
                # interaction with `other < spec`
                for other in 1:spec-1
                    local Wprime = Wprimes[spec][other]
                    local w::eltype(dx.x[1]) = 0
                    for j in 1:length(x.x[other])
                        w += Wprime(x.x[other][j] - x.x[spec][i])
                    end
                    v += w / (length(x.x[other]) - 1)
                end
                # interaction with `other > spec`
                for other in spec+1:N
                    local Wprime = Wprimes[spec][other]
                    local w::eltype(dx.x[1]) = 0
                    for j in 1:length(x.x[other])
                        w += Wprime(x.x[other][j] - x.x[spec][i])
                    end
                    v += w / (length(x.x[other]) - 1)
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

export sampled_interaction, integrated_interaction

function sampled_interaction(x::Real, Wprime, ys::AbstractVector{<:Real})
    sum(Wprime(y - x) for y in ys) / (length(ys) - 1)
end

function sampled_interaction(t::Real, x::Real, Wprime, ys::AbstractVector{<:Real})
    sum(Wprime(t, y - x) for y in ys) / (length(ys) - 1)
end

function sampled_interaction(x::Real; Wprime, particles::AbstractVector{<:Real})
    sum(Wprime(p - x) for p in particles) / (length(particles) - 1)
end

function sampled_interaction(t::Real, x::Real; Wprime, particles::AbstractVector{<:Real})
    sum(Wprime(t, p - x) for p in particles) / (length(particles) - 1)
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


function integrated_interaction_(x::Real, W, ys::AbstractVector{<:Real}, dens::AbstractVector{<:Real})
    sum(i -> (dens[i] - dens[i+1]) * W(ys[i] - x), eachindex(ys))
end

function integrated_interaction_(t::Real, x::Real, W, s::Integer, dens::AbstractArray{<:Real,3}, ys::AbstractVector{<:Real})
    sum(i -> (dens[s, 1, i] - dens[s, 2, i]) * W(t, ys[i] - x), eachindex(ys))
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
                    v += sampled_interaction(Wprimes[spec][other], x.x[other], x.x[spec][i])
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

export make_velocities_time
function make_velocities_time(
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
                v::F = Vs[spec](t, x.x[spec][i])
                for other in 1:N
                    v += sampled_interaction(t, Wprimes[spec][other], x.x[other], x.x[spec][i])
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


export SampledInteraction
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


export IntegratedInteraction
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


export velocities, velocities!

function velocities(
    x::ArrayPartition{F, T},
    p::Union{SampledInteraction, IntegratedInteraction},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}}
}
    dx = similar(x)
    velocities!(dx, x, p, t)
    dx
end

function velocities!(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::SampledInteraction{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
    dens = pwc_densities(x.x...)
    for spec in 1:N
        dx.x[spec] .= p.Vs[spec].(t, x.x[spec])
        for other in 1:N, i in eachindex(x.x[spec])
            dx.x[spec][i] += sampled_interaction(t, x.x[spec][i], p.Wprimes[spec][other], x.x[other])
        end
        # for other in 1:N
        #     dx.x[spec] .+= sampled_interaction.(t, x.x[spec]; Wprime=p.Wprimes[spec][other], particles=x.x[other])
        # end
        d = dens[spec]
        for i in eachindex(x.x[spec])
            if dx.x[spec][i] < 0
                mob = p.mobilities[spec](d[:, 1, i]...)
            else
                mob = p.mobilities[spec](d[:, 2, i]...)
            end
            dx.x[spec][i] *= mob
        end
    end
end

function velocities!(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::IntegratedInteraction{N, TVs, TWs, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWs, Tmobilities
}
    dens = pwc_densities(x.x...)
    dens_diff = similar(x)
    for s in 1:N
        dens_diff.x[s] .= dens[s][s, 1, :] .- dens[s][s, 2, :]
    end
    for spec in 1:N
        dx.x[spec] .= p.Vs[spec].(t, x.x[spec])
        for other in 1:N, i in eachindex(x.x[spec])
            dx.x[spec][i] += integrated_interaction(t, x.x[spec][i], p.Ws[spec][other], x.x[other], dens_diff.x[other])
        end
        d = dens[spec]
        for i in eachindex(x.x[spec])
            if dx.x[spec][i] < 0
                mob = p.mobilities[spec](d[:, 1, i]...)
            else
                mob = p.mobilities[spec](d[:, 2, i]...)
            end
            dx.x[spec][i] *= mob
        end
    end
end


export velocities_gen, velocities_gen!

function velocities_gen(
    x::ArrayPartition{F, T},
    p::Union{SampledInteraction, IntegratedInteraction},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}}
}
    dx = similar(x)
    velocities_gen!(dx, x, p, t)
    dx
end

@generated function velocities_gen!(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::SampledInteraction{N, TVs, TWprimes, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWprimes, Tmobilities
}
quote
    dens = pwc_densities(x.x...)
    $((quote
        dx.x[$spec] .= p.Vs[$spec].(t, x.x[$spec])
        $((quote
            for i in eachindex(x.x[$spec])
                dx.x[$spec][i] += sampled_interaction(t, x.x[$spec][i], p.Wprimes[$spec][$other], x.x[$other])
            end
        end for other in 1:N)...)
        d = dens[$spec]
        for i in eachindex(x.x[$spec])
            if dx.x[$spec][i] < 0
                mob = p.mobilities[$spec]($((:(d[$j, 1, i]) for j in 1:N)...))
            else
                mob = p.mobilities[$spec]($((:(d[$j, 2, i]) for j in 1:N)...))
            end
            dx.x[$spec][i] *= mob
        end
    end for spec in 1:N)...)
end
end

@generated function velocities_gen!(
    dx::ArrayPartition{F, T},
    x::ArrayPartition{F, T},
    p::IntegratedInteraction{N, TVs, TWs, Tmobilities},
    t
) where {
    F,
    T <: Tuple{Vararg{AbstractVector{<:Real}}},
    N, TVs, TWs, Tmobilities
}
quote
    dens = pwc_densities(x.x...)
    dens_diff = similar(x)
    for s in 1:N
        dens_diff.x[s] .= dens[s][s, 1, :] .- dens[s][s, 2, :]
    end
    # $((:(dens_diff.x[$s] .= dens[$s][$s, 1, :] .- dens[$s][$s, 2, :]) for s in 1:N)...)
    $((quote
        dx.x[$spec] .= p.Vs[$spec].(t, x.x[$spec])
        $((quote
            for i in eachindex(x.x[$spec])
                dx.x[$spec][i] += integrated_interaction(t, x.x[$spec][i], p.Ws[$spec][$other], x.x[$other], dens_diff.x[$other])
            end
        end for other in 1:N)...)
        d = dens[$spec]
        for i in eachindex(x.x[$spec])
            if dx.x[$spec][i] < 0
                mob = p.mobilities[$spec]($((:(d[$j, 1, i]) for j in 1:N)...))
            else
                mob = p.mobilities[$spec]($((:(d[$j, 2, i]) for j in 1:N)...))
            end
            dx.x[$spec][i] *= mob
        end
    end for spec in 1:N)...)
end
end
