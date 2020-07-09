module Particles

using RecursiveArrayTools
using DocStringExtensions

include("utils.jl")
include("densities.jl")
include("velocities.jl")
include("plotting.jl")
include("examples.jl")

export empty_like
export pwc_density, pwc_densities
export make_velocity, make_velocities
export plot_density, plot_density!, plot_linear_density!

end # module
