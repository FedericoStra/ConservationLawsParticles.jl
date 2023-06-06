export PiecewiseConstantFunction
export lpnorm

struct PiecewiseConstantFunction{X,Y}
    cuts::Vector{X}
    values::Vector{Y}
end

Base.map(f, pcf::PiecewiseConstantFunction) = PiecewiseConstantFunction(pcf.cuts, map(f, pcf.values))

Base.abs(pcf::PiecewiseConstantFunction) = map(abs, pcf)

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

function lpnorm(f::PiecewiseConstantFunction, p)
    if isone(p)
        sum(i -> abs(f.values[i]) * (f.cuts[i]-f.cuts[i-1]), 2:length(f.values)-1)
    elseif isinf(p)
        maximum(abs, f.values)
    else
        sum(i -> abs(f.values[i])^p * (f.cuts[i]-f.cuts[i-1]), 2:length(f.values)-1)^(1/p)
    end
end
