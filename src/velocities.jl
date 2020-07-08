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
                v += Vs[spec](x[i])
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

function total_interaction_(Wprime, ys::AbstractVector{<:Real}, x::Real)
    T = promote_type(eltype(ys), typeof(x))
    w::T = 0
    for y in ys
        w += Wprime(y - x)
    end
    w / length(ys)
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
