"""
$(TYPEDSIGNATURES)

    pwc_density(x::AbstractVector)

Returns the piecewise constant probability density from the quantile particle positions.

Let the quantile particles be at positions `x₀, x₁, …, xₙ` (remember that in Julia they are actually `x[1], x[2], ..., x[n], x[n+1]`).
The returned array `R` is such that `R[1] = R[n+2] = 0` and `R[i] = 1 / (N * (x[i] - x[i-1]))` for the intermediate indices (this formula holds also for the first and last entry if we assume `x[0] = -∞` and `x[n+2] = ∞`).

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
function pwc_density(x::AbstractVector)
    len = length(x)
    R = Array{float(eltype(x))}(undef, len+1)
    R[1] = R[end] = 0
    d = 1 / (len-1)
    for i in 2:len
        R[i] = d / (x[i] - x[i-1])
    end
    R
end
