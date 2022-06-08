using Distributed
using Graphs
using JuMP, GLPK
using Graphs.Parallel

∇(f, e::Te where Te <: AbstractEdge) = f[e.dst] - f[e.src]

function κ(
    G::Tg where Tg <: AbstractGraph, 
    e::Te where Te <: AbstractEdge,
    D::AbstractArray,
    Δ::AbstractArray; 
    optimizer = GLPK.Optimizer
    )

    here = unique( [neighbors(G, e.src)..., neighbors(G, e.dst)...] )

    d = D[here, here]
    δ = Δ[here, here]
    

    model = Model(optimizer)
    @variable(model, f[here])
    @constraint(model, .- d .≤ [f[x] - f[y] for x ∈ here, y ∈ here] .≤ d)
    @constraint(model, ∇(f, e) == 1)
    @objective(model, Min, ∇(δ * f, e))

    optimize!(model)

    return objective_value(model)

end

function κ(
    G::AbstractGraph,
    E::Union{Graphs.AbstractEdgeIter, Vector{T} where T <: AbstractEdge} = edges(G);
    optimizer = GLPK.Optimizer,
    parallel = false
    )

    Δ = laplacian_matrix(G)
    D = (parallel ? Parallel : Graphs).floyd_warshall_shortest_paths(squash(G)).dists

    return (parallel ? pmap : map)(e -> κ(G, e, D, Δ; optimizer = optimizer), E)

end