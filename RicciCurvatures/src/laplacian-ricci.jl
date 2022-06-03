using Distributed
using Graphs
using JuMP, GLPK

∇(f, e::Te where Te <: AbstractEdge) = f[e.dst] - f[e.src]

function κ(
    G::Tg where Tg <: AbstractGraph, 
    e::Te where Te <: AbstractEdge; 
    optimizer = GLPK.Optimizer
    )

    here = unique( [neighbors(G, e.src)..., neighbors(G, e.dst)...] )
    n = length(here)
    g, _ = induced_subgraph(G, here)

    Δ = laplacian_matrix(g)
    D = floyd_warshall_shortest_paths(squash(g)).dists

    model = Model(optimizer)
    @variable(model, f[1:n])
    @constraint(model, .- D .≤ [f[x] - f[y] for x ∈ 1:n, y ∈ 1:n] .≤ D)
    @constraint(model, ∇(f, e) == 1)
    @objective(model, Min, ∇(Δ * f, e))

    optimize!(model)

    return objective_value(model)

end

function κ(
    G::AbstractGraph,
    E::Union{Graphs.AbstractEdgeIter, Vector{T} where T <: AbstractEdge} = edges(G);
    optimizer = GLPK.Optimizer,
    parallel = false
    )

    return (parallel ? pmap : map)(e -> κ(G, e; optimizer = optimizer), E)

end