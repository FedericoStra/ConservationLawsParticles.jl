export pwc_density, pwc_densities, pwc_densities!

"""
    pwc_density(x::AbstractVector)

Returns the piecewise constant probability density from the quantile particle positions.

Let the quantile particles be at positions `x₀, x₁, …, xₙ` (remember that in Julia
they are actually `x[1], x[2], ..., x[n], x[n+1]`).
The returned array `R` is such that `R[1] = R[n+2] = 0` and `R[i] = 1 / (n * (x[i] - x[i-1]))`
for the intermediate indices (this formula holds also for the first and last entry
if we assume `x[0] = -∞` and `x[n+2] = ∞`).

# Examples

```jldoctest
julia> pwc_density([0, 1, 3])
4-element Vector{Float64}:
 0.0
 0.5
 0.25
 0.0
```
"""
function pwc_density(x::AbstractVector{<:Real})
    len = length(x)
    R = Array{float(eltype(x))}(undef, len+1)
    R[1] = R[end] = 0
    d = 1 / (len-1)
    for i in 2:len
        R[i] = d / (x[i] - x[i-1])
    end
    R
end

"""
    pwc_density(x::AbstractVector, y::AbstractVector)

The returned `x_dens` is indexed as `x_dens[x_or_y, left_or_right, i]`.
"""
function pwc_density(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    x_len, y_len = length(x), length(y)
    T = promote_type(float(eltype(x)), float(eltype(y)))
    x_dens = Array{T, 3}(undef, (2, 2, x_len))
    y_dens = Array{T, 3}(undef, (2, 2, y_len))
    # indices of current particles being examined
    x_i::Int, y_i::Int = 1, 1
    # current densities
    x_d::T, y_d::T = zero(T), zero(T)
    # normalization factors
    x_f::T, y_f::T = 1 / (x_len - 1), 1 / (y_len - 1)
    # @debug "Entering loop" x_len y_len
    @inbounds while x_i <= x_len || y_i <= y_len
        if y_i > y_len || x_i <= x_len && x[x_i] < y[y_i]
            # @debug "Iteration x < y" x_i y_i
            x_dens[1, 1, x_i] = x_d
            x_dens[2, 1, x_i] = y_d
            if x_i < x_len
                x_d = x_f / (x[x_i + 1] - x[x_i])
            else
                x_d = zero(T)
            end
            x_dens[1, 2, x_i] = x_d
            x_dens[2, 2, x_i] = y_d
            x_i += 1
        elseif x_i > x_len || y_i <= y_len && y[y_i] < x[x_i]
            # @debug "Iteration y < x" x_i y_i
            y_dens[1, 1, y_i] = x_d
            y_dens[2, 1, y_i] = y_d
            if y_i < y_len
                y_d = y_f / (y[y_i + 1] - y[y_i])
            else
                y_d = zero(T)
            end
            y_dens[1, 2, y_i] = x_d
            y_dens[2, 2, y_i] = y_d
            y_i += 1
        else
            # @debug "Iteration x == y" x_i y_i
            x_dens[1, 1, x_i] = x_d
            x_dens[2, 1, x_i] = y_d
            y_dens[1, 1, y_i] = x_d
            y_dens[2, 1, y_i] = y_d
            if x_i < x_len
                x_d = x_f / (x[x_i + 1] - x[x_i])
            else
                x_d = zero(T)
            end
            if y_i < y_len
                y_d = y_f / (y[y_i + 1] - y[y_i])
            else
                y_d = zero(T)
            end
            x_dens[1, 2, x_i] = x_d
            x_dens[2, 2, x_i] = y_d
            y_dens[1, 2, y_i] = x_d
            y_dens[2, 2, y_i] = y_d
            x_i += 1
            y_i += 1
        end
    end
    (x_dens, y_dens)
end


"""
    pwc_densities(xs::AbstractVector{<:Real}...)
    pwc_densities(xs::Tuple{Vararg{AbstractVector{<:Real}}})

Given a list of `n` families of particles of lengths `l₁, …, lₙ`,
returns an `n`-tuple whose `s`-th entry is an `n×2×lₛ` array `densₛ`.

This array contains the information about the densities of all the species
around the particles of the `s`-th species.
This array is indexed as `[o,side,i]`, where `o∈{1, …, n}` is the index of the other
species whose density we want to know, `side` is either `1` if we want to know
the left density or `2` for the right density, and `i∈{1, …, lₛ}` is the index
of the particle.

In other words, `densₛ[o,side,i]` is the density of the `o`-th species on the left (`side=1`)
or right (`side=2`) of the `i`-th particle of the `s`-th species.

See also [`pwc_densities!`](@ref) for an inplace version.

# Examples

```jldoctest
julia> pwc_densities([0,1,2], [0,2,3])
([0.0 0.5; 0.0 0.25;;; 0.5 0.5; 0.25 0.25;;; 0.5 0.0; 0.25 0.5], [0.0 0.5; 0.0 0.25;;; 0.5 0.0; 0.25 0.5;;; 0.0 0.0; 0.5 0.0])
```
"""
function pwc_densities end

pwc_densities(xs::AbstractVector{<:Real}...) = pwc_densities(xs)

function pwc_densities(xs::Tuple{Vararg{AbstractVector{<:Real}}})
    S = nfields(xs)
    T = promote_type(float.(eltype.(xs))...)
    len = length.(xs)
    dens = map(x -> new_undef_densities(T, S, length(x)), xs)::NTuple{S, Array{T, 3}}
    pwc_densities!(dens, xs)
end


"""
    pwc_densities!(dens, xs::AbstractVector{<:Real}...)
    pwc_densities!(dens, xs::Tuple{Vararg{AbstractVector{<:Real}}})

This is an inplace version of [`pwc_densities`](@ref) that writes the result in `dens`.

See the documentation of [`pwc_densities`](@ref) for an explanation.
"""
function pwc_densities! end

pwc_densities!(dens, xs::AbstractVector{<:Real}...) = pwc_densities!(dens, xs)

function pwc_densities!(dens, xs::Tuple{Vararg{AbstractVector{<:Real}}})
    S = nfields(xs)
    T = promote_type(float.(eltype.(xs))...)
    len = length.(xs)
    # indices of current particles being examined
    ind = ones(Int, S)
    # current densities
    ds = zeros(T, S)
    # normalization factors
    fs::NTuple{S, T} = 1 ./ (len .- 1)
    # @debug "Entering loop" T len ds fs
    i_min = Int[0]
    @inbounds while true
        n_mins::Int = 0
        for i in 1:S
            if ind[i] <= len[i]
                if n_mins == 0 || xs[i][ind[i]] < xs[i_min[1]][ind[i_min[1]]]
                    i_min[1] = i
                    n_mins = 1
                elseif xs[i][ind[i]] <= xs[i_min[1]][ind[i_min[1]]]
                    n_mins += 1
                    if n_mins <= length(i_min)
                        i_min[n_mins] = i
                    else
                        push!(i_min, i)
                    end
                end
            end
        end
        if n_mins == 0
            break
        end
        i_mins = @view i_min[1:n_mins]
        # let ind = (ind...,), i_min = (i_min...,), i_mins = (i_mins...,)
        #     @debug "Iteration" ind i_min n_mins i_mins
        # end
        for i in i_mins
            # @debug "Assigning left" i ind[i] ds
            # dens[i][:, 1, ind[i]] = ds
            for s in 1:S
                dens[i][s, 1, ind[i]] = ds[s]
            end
        end
        for i in i_mins
            if ind[i] < length(xs[i])
                ds[i] = fs[i] / (xs[i][ind[i]+1] - xs[i][ind[i]])
            else
                ds[i] = zero(T)
            end
            # ds[i] = right_local_density(xs[i], ind[i], fs[i])
        end
        for i in i_mins
            # @debug "Assigning right" i ind[i] ds
            # dens[i][:, 2, ind[i]] = ds
            for s in 1:S
                dens[i][s, 2, ind[i]] = ds[s]
            end
        end
        # ind[i_mins] .+= 1
        for i in i_mins
            ind[i] += 1
        end
    end
    dens
end

new_undef_densities(T::Type, S::Int, len::Int) = Array{T, 3}(undef, (S, 2, len))
