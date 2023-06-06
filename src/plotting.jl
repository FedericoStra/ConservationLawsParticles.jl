export densityplot, densityplot!
export plot_density, plot_density!, plot_linear_density!

export pcfplot, pcfplot!
export plot_pcf, plot_pcf!

using RecipesBase


@userplot DensityPlot

@recipe function _(dp::DensityPlot)
    [(
        X = Array{eltype(x)}(undef, 2*length(x));
        Y = Array{Float64}(undef, 2*length(x));
        R = pwc_density(x);
        for i in eachindex(x)
            X[2i-1] = X[2i] = x[i]
            Y[2i-1] = R[i]
            Y[2i]   = R[i+1]
        end;
        (X, Y)
    ) for x in dp.args]
end

plot_density(args...; kw...) = densityplot(args...; kw...)
plot_density!(args...; kw...) = densityplot!(args...; kw...)

function plot_linear_density!(p, x; kwargs...)
    X = (x[1:end-1] + x[2:end]) / 2
    Y = pwc_density(x)[2:end-1]
    plot!(p, X, Y; kwargs...)
end


@userplot PcfPlot

@recipe function _(pcfp::PcfPlot)
    [(
        X = repeat(f.cuts, inner=2);
        Y = repeat(f.values, inner=2)[2:end-1];
        X[1] = f.cuts[1];
        (X, Y)
    ) for f in pcfp.args]
end

plot_pcf(args...; kwargs...) = pcfplot(args...; kwargs...)
plot_pcf!(args...; kwargs...) = pcfplot!(args...; kwargs...)
