export Diffusion, MinDiffusion, MeanDiffusion, SimpleDiffusion, diffuse!

struct Diffusion{D}
    diffusion::D
end

struct MinDiffusion{D}
    diffusion::D
end

struct MeanDiffusion{D}
    diffusion::D
end

struct SimpleDiffusion{D}
    diffusion::D
end

function diffuse!(dx, x, dens, diffusion::Nothing) end

function diffuse!(dx, x, dens, diffusion::Real)
    diffuse!(dx, x, dens, Diffusion(diffusion))
end

function diffuse!(dx, x, dens, diffusion::Diffusion{<:Real})
    diffuse!(dx, x, dens, MinDiffusion(diffusion.diffusion))
end

function diffuse!(dx, x, dens, diffusion::MinDiffusion{<:Real})
    for i in eachindex(dx)
        δdens = dens[i+1] - dens[i]
        if i == 1
            δx = x[2] - x[1]
            ρ = dens[2]
        elseif i == length(dx)
            δx = x[end] - x[end-1]
            ρ = dens[end-1]
        else
            δx = (x[i+1] - x[i-1]) / 2
            ρ = min(dens[i], dens[i+1])
        end
        dx[i] -= diffusion.diffusion * δdens / δx / ρ
    end
end

function diffuse!(dx, x, dens, diffusion::MeanDiffusion{<:Real})
    for i in eachindex(dx)
        δdens = dens[i+1] - dens[i]
        if i == 1
            δx = x[2] - x[1]
        elseif i == length(dx)
            δx = x[end] - x[end-1]
        else
            δx = (x[i+1] - x[i-1]) / 2
        end
        ρ = (dens[i] + dens[i+1]) / 2
        dx[i] -= diffusion.diffusion * δdens / δx / ρ
    end
end

function diffuse!(dx, x, dens, diffusion::SimpleDiffusion{<:Real})
    for i in eachindex(dx)
        δdens = dens[i+1] - dens[i]
        dx[i] -= diffusion.diffusion * (length(dx) - 1) * δdens
    end
end
