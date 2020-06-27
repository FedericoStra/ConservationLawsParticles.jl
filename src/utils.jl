"""
    empty_like(x::AbstractArray)

Creates an empty array with the same `size` and `eltype` of `x`.

# Examples
```julia-repl
julia> empty_like([1 2 3; 4 5 6])
2×3 Array{Int64,2}:
 0  140456850103040                1
 1                3  140456850103041
```
"""
empty_like(x::AbstractArray) = Array{eltype(x)}(undef, size(x))
