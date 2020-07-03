"""
$(TYPEDSIGNATURES)

    pwc_density(x::AbstractVector)

Returns the piecewise constant probability density from the quantile particle positions.

Let the quantile particles be at positions `x₀, x₁, …, xₙ` (remember that in Julia
they are actually `x[1], x[2], ..., x[n], x[n+1]`).
The returned array `R` is such that `R[1] = R[n+2] = 0` and `R[i] = 1 / (N * (x[i] - x[i-1]))`
for the intermediate indices (this formula holds also for the first and last entry
if we assume `x[0] = -∞` and `x[n+2] = ∞`).

# Examples

```jldoctest; setup = :(using Particles)
julia> pwc_density([0, 1, 3])
4-element Array{Float64,1}:
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
    while x_i <= x_len || y_i <= y_len
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
    pwc_densities(xs::AbstractVector{<:AbstractVector})
"""
# function pwc_densities(xs::AbstractVector{<:AbstractVector})
#     len = length.(xs)
#     T = float(promote_type(eltype.(xs)...))
# end

"""
    pwc_densities(xs::AbstractVector...)
"""
function pwc_densities(xs::Vararg{AbstractVector{<:Real}, N}) where N
    N::Int
    T = promote_type(float.(eltype.(xs))...)
    len = length.(xs)
    dens = map(x -> new_undef_densities(T, N, length(x)), xs)::NTuple{N, Array{T, 3}}
    ind = ones(Int, N)
    ds = zeros(T, N)
    fs::NTuple{N, T} = 1 ./ (len .- 1)
    # @debug "Entering loop" T len ds fs
    while true
        i_min = Int[]
        for i in 1:N
            if ind[i] <= len[i]
                if isempty(i_min) || xs[i][ind[i]] < xs[i_min[1]][ind[i_min[1]]]
                    i_min = [i]
                elseif xs[i][ind[i]] <= xs[i_min[1]][ind[i_min[1]]]
                    push!(i_min, i)
                end
            end
        end
        if isempty(i_min)
            break
        end
        # let ind = (ind...,), i_min = (i_min...,)
        #     @debug "Iteration" ind i_min
        # end
        for i in i_min
            # @debug "Assigning left" i ind[i] ds
            dens[i][:, 1, ind[i]] .= ds
        end
        for i in i_min
            if ind[i] < length(xs[i])
                ds[i] = fs[i] / (xs[i][ind[i]+1] - xs[i][ind[i]])
            else
                ds[i] = zero(T)
            end
            # ds[i] = right_local_density(xs[i], ind[i], fs[i])
        end
        for i in i_min
            # @debug "Assigning right" i ind[i] ds
            dens[i][:, 2, ind[i]] .= ds
        end
        ind[i_min] .+= 1
    end
    dens
end

new_undef_densities(T::Type, N::Int, len::Int) = Array{T, 3}(undef, (N, 2, len))

function right_local_density(x::AbstractVector, i::Integer, f::Number)
    if i < length(x)
        f / (x[i+1] - x[i])
    else
        0
    end
end
