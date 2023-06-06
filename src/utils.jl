using SpecialFunctions

export gaussian_particles, @time_independent


"""
    gaussian_particles(width::Real, number::Integer)

Generates `number` particles defining a Gaussian distribution with variance `1/2` truncated
to the interval `[-width, width]`.
"""
function gaussian_particles(width::Real, number::Integer)
    width > 0 && number > 1 || error("requirements: width > 0 and number > 1")
    e = erf(width)
    erfinv.(range(-e, e, length=number))
end


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
