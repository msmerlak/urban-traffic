using Test
using RicciCurvatures
using Graphs

@testset "Square" begin
    V = 4
    E = 4
    G = SimpleGraph(V);
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 3)
    add_edge!(G, 3, 4)
    add_edge!(G, 4, 1)

    @test κ(G) ≈ fill(V/2, E)
end

@testset "Complete graphs" begin
    for V in 2:10
        E = Int(V * (V-1) / 2)
        G = complete_graph(V);
        @test κ(G) ≈ fill(V, E)
    end
end

