export PiecewiseConstantFunction
export plot_pcf, plot_pcf!

struct PiecewiseConstantFunction{X,Y}
    cuts::Vector{X}
    values::Vector{Y}
end

Base.broadcast(f, pcf::PiecewiseConstantFunction) = PiecewiseConstantFunction(pcf.cuts, broadcast(f, pcf.values))

Base.map(f, pcf::PiecewiseConstantFunction) = PiecewiseConstantFunction(pcf.cuts, map(f, pcf.values))

function Base.broadcast(op, f::PiecewiseConstantFunction, g::PiecewiseConstantFunction)
    cuts = promote_type(eltype(f.cuts), eltype(g.cuts))[]
    values = promote_type(eltype(f.values), eltype(g.values))[]
    i, j = 1, 1
    len_f, len_g = length(f.cuts), length(g.cuts)
    while i <= len_f || j <= len_g
        if j > len_g || i <= len_f && f.cuts[i] < g.cuts[j]
            push!(cuts, f.cuts[i])
            push!(values, op(f.values[i], g.values[j]))
            i += 1
        elseif i > len_f || j <= len_g && g.cuts[j] < f.cuts[i]
            push!(cuts, g.cuts[j])
            push!(values, op(f.values[i], g.values[j]))
            j += 1
        else
            push!(cuts, f.cuts[i])
            push!(values, op(f.values[i], g.values[j]))
            i += 1
            j += 1
        end
    end
    push!(values, op(f.values[i], g.values[j]))
    PiecewiseConstantFunction(cuts, values)
end

for op in [:+, :-, :*, :/, :^]
@eval function Base.$op(f::PiecewiseConstantFunction, g::PiecewiseConstantFunction)
    cuts = promote_type(eltype(f.cuts), eltype(g.cuts))[]
    values = promote_type(eltype(f.values), eltype(g.values))[]
    i, j = 1, 1
    len_f, len_g = length(f.cuts), length(g.cuts)
    while i <= len_f || j <= len_g
        if j > len_g || i <= len_f && f.cuts[i] < g.cuts[j]
            push!(cuts, f.cuts[i])
            push!(values, $op(f.values[i], g.values[j]))
            i += 1
        elseif i > len_f || j <= len_g && g.cuts[j] < f.cuts[i]
            push!(cuts, g.cuts[j])
            push!(values, $op(f.values[i], g.values[j]))
            j += 1
        else
            push!(cuts, f.cuts[i])
            push!(values, $op(f.values[i], g.values[j]))
            i += 1
            j += 1
        end
    end
    push!(values, $op(f.values[i], g.values[j]))
    PiecewiseConstantFunction(cuts, values)
end
end

function plot_pcf(f::PiecewiseConstantFunction; kwargs...)
    x = f.cuts
    R = f.values
    X = Array{eltype(x)}(undef, 2*length(x)-2)
    Y = Array{Float64}(undef, 2*length(x)-2)
    X[1] = x[1]
    for i in 2:length(x)-1
        X[2i-2] = X[2i-1] = x[i]
        Y[2i-3] = Y[2i-2] = R[i]
    end
    X[end] = x[end]
    Y[end-1] = Y[end] = R[end-1]
    plot(X, Y; kwargs...)
end
    
function plot_pcf!(p, f::PiecewiseConstantFunction; kwargs...)
    x = f.cuts
    R = f.values
    X = Array{eltype(x)}(undef, 2*length(x)-2)
    Y = Array{Float64}(undef, 2*length(x)-2)
    X[1] = x[1]
    for i in 2:length(x)-1
        X[2i-2] = X[2i-1] = x[i]
        Y[2i-3] = Y[2i-2] = R[i]
    end
    X[end] = x[end]
    Y[end-1] = Y[end] = R[end-1]
    plot!(p, X, Y; kwargs...)
end
