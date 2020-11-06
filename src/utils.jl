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
function empty_like end
@deprecate empty_like(x::AbstractArray) similar(x)
# empty_like(x::AbstractArray) = Array{eltype(x)}(undef, size(x))

using SpecialFunctions
export gaussian_particles
function gaussian_particles(width::Real, number::Integer)
    width > 0 && number > 1 || error("requirements: width > 0 and number > 1")
    e = erf(width)
    erfinv.(range(-e, e, length=number))
end


export @time_independent

"""
Automatically define a method which takes time as first argument and discards it.

# Examples
The definition
```julia
@time_independent V(x) = -x^3
```
is equivalent to
```julia
V(x) = -x^3
V(t, x) = V(x)
```

This also works with more than one argument, for instance
```julia
@time_independent V(x₁, x₂) = x₁ * x₂
```
is equivalent to
```julia
V(x₁, x₂) = x₁ * x₂
V(t, x₁, x₂) = V(x₁, x₂)
```
"""
macro time_independent(code)
    @assert code isa Expr
    @assert code.head == :(=) || code.head == :function
    lhs = code.args[1]
    @assert lhs.head == :call
    f = lhs.args[1]
    quote
        $(esc(code))
        $(esc(f))(t, args...) = $(esc(f))(args...)
    end
end
