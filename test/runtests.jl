using SafeTestsets
using Test
using TestSetExtensions

@testset ExtendedTestSet "All tests" begin

    @safetestset "Aqua tests" include("Aqua.jl")

    @testset "ConservationLawsParticles.jl" begin
        @safetestset "Densities"  include("test_densities.jl")
        @safetestset "Velocities" include("test_velocities.jl")
        @safetestset "Utils"      include("test_utils.jl")
    end

end
