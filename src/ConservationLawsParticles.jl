module ConservationLawsParticles

using RecursiveArrayTools
using DocStringExtensions

include("utils.jl")
include("densities.jl")
include("models.jl")
include("interactions.jl")
include("diffusion.jl")
include("velocities.jl")
include("plotting.jl")
include("examples.jl")

end # module
