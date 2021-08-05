using Plots

function plot_density!(p, x; kwargs...)
    X = Array{eltype(x)}(undef, 2*length(x))
    Y = Array{Float64}(undef, 2*length(x))
    R = pwc_density(x)
    for i in eachindex(x)
        X[2i-1] = X[2i] = x[i]
        Y[2i-1] = R[i]
        Y[2i]   = R[i+1]
    end
    plot!(p, X, Y; kwargs...)
end

function plot_density(x; kwargs...)
    p = plot(; kwargs...)
    plot_density!(p, x; kwargs...)
end

function plot_linear_density!(p, x; kwargs...)
    X = (x[1:end-1] + x[2:end]) / 2
    Y = pwc_density(x)[2:end-1]
    plot!(p, X, Y; kwargs...)
end
