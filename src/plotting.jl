using Plots

function plot_density!(p, x; kwargs...)
    X = Array{eltype(x)}(undef, 2*length(x)-2)
    Y = Array{Float64}(undef, 2*length(x)-2)
    R = pwc_density(x)
    X[1] = x[1]
    for i in 2:length(x)-1
        X[2i-2] = X[2i-1] = x[i]
        Y[2i-3] = Y[2i-2] = R[i]
    end
    X[end] = x[end]
    Y[end-1] = Y[end] = R[end-1]
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
