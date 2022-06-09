using Test
using RicciCurvatures
using Graphs

@testset "Hypercube graphs" begin
    K₂ = complete_graph(2)
    H = copy(K₂)
    for n in 2:6
        H = cartesian_product(H, K₂)
        @test κ(H) ≈ fill(2/n, length(edges(H)))
    end
end

@testset "Cycle graphs" begin
    @test κ(cycle_graph(3)) ≈ fill(3/2, 3)
    @test κ(cycle_graph(4)) ≈ fill(1, 4)
    @test κ(cycle_graph(5)) ≈ fill(1/2, 5)
    for V in 6:10
        E = V
        @test κ(cycle_graph(V)) ≈ zeros(E)
    end
end


@testset "Complete graphs" begin
    for V in 2:10
        G = complete_graph(V)
        E = Int(V*(V-1)/2)
        @test κ(G) ≈ fill(V/(V-1), E)
    end
end

