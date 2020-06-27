module Particles

using DocStringExtensions

include("densities.jl")
include("utils.jl")
include("plotting.jl")

export empty_like, pwc_density
export plot_density, plot_density!, plot_linear_density!

end # module
