using ConservationLawsParticles
using RecursiveArrayTools

@testset "pwc_density" begin
    # test a few specific examples
    @testset "examples" begin
        @test pwc_density([0.]) == [0.0, 0.0]
        @test pwc_density([0., 1.]) == [0.0, 1.0, 0.0]
        @test pwc_density([0, 1, 3]) == [0.0, 0.5, 0.25, 0.0]
    end

    # test 1-S pwc_density for sum == 1
    @testset "1-S sum == 1" begin
        @testset for len in [2, 5, 10, 50, 100, 500, 1000]
            x = sort(randn(len))
            d = pwc_density(x)
            @test sum(d[2:end-1] .* diff(x)) ≈ 1.  atol=100eps()
        end
    end

    # test 2-S pwc_density for sum == 1
    @testset "2-S sum == 1" begin
        @testset for len in [2, 5, 10, 50, 100, 500, 1000]
            x = ArrayPartition((sort(randn(len)) for _ in 1:2)...)
            d = pwc_density(x.x...)
            for s in 1:2
                @test sum(d[s][s, 1, 2:end] .* diff(x.x[s])) ≈ 1.  atol=100eps()
            end
        end
    end

    # test left/right consistency of `pwc_density` with 2 species
    @testset "left/right" begin
        @testset for len in [2, 3, 10, 100]
            x = ArrayPartition((sort(randn(len)) for _ in 1:2)...)
            d = pwc_density(x.x...)
            for s in 1:2
                @test d[s][s, 1, 2:end] == d[s][s, 2, 1:end-1]
            end
        end
    end

    # test consistency of `pwc_density` with 1 and 2 species
    @testset "1-S/2-S" begin
        @testset for len in [2, 3, 10, 100]
            x = ArrayPartition((sort(randn(len)) for _ in 1:2)...)
            d = pwc_density(x.x...)
            for s in 1:2
                @test d[s][s, 1, 2:end] == pwc_density(x.x[s])[2:end-1]
            end
        end
    end
end

@testset "pwc_densities" begin
    # test left/right consistency of `pwc_densities` with N species
    @testset "left/right" begin
        @testset for N in 1:5
            @testset for len in [2, 3, 10, 100]
                x = ArrayPartition((sort(randn(len)) for _ in 1:N)...)
                d = pwc_densities(x.x...)
                for s in 1:N
                    @test d[s][s, 1, 2:end] == d[s][s, 2, 1:end-1]
                end
            end
        end
    end

    # test N-S pwc_densities for sum == 1
    @testset "sum == 1" begin
        @testset for N in 1:5
            @testset for len in [2, 5, 10, 50, 100, 500, 1000]
                x = ArrayPartition((sort(randn(len)) for _ in 1:N)...)
                d = pwc_densities(x.x...)
                for s in 1:N
                    @test sum(d[s][s, 1, 2:end] .* diff(x.x[s])) ≈ 1.  atol=100eps()
                end
            end
        end
    end

    # test consistency of `pwc_densities` with 1 and N species
    @testset "1-S/N-S" begin
        @testset for N in 1:5
            @testset for len in [2, 3, 10, 100]
                x = ArrayPartition((sort(randn(len)) for _ in 1:N)...)
                d = pwc_densities(x.x...)
                for s in 1:N
                    @test d[s][s:s, :, :] == pwc_densities(x.x[s])[1]
                end
            end
        end
    end
end

# test 1-S/2-S pwc_density == pwc_densities
@testset "pwc_density == pwc_densities" begin
    @testset "1-S" begin
        @testset for len in [2, 3, 10, 100]
            x = sort(randn(len))
            @test pwc_density(x)[2:end-1] == pwc_densities(x)[1][1, 1, 2:end]
        end
    end

    @testset "2-S" begin
        @testset for len in [2, 3, 10, 100]
            x = ArrayPartition((sort(randn(len)) for _ in 1:2)...)
            @test pwc_density(x.x...) == pwc_densities(x.x...)
        end
    end

    @testset "2-S (common points)" begin
        common = collect(-1:0.1:1)
        @testset for lengths in [(2,3), (5,10), (50,100)]
            x = ArrayPartition((sort(vcat(randn(l), common)) for l in lengths)...)
            @test pwc_density(x.x...) == pwc_densities(x.x...)
        end
    end
end
