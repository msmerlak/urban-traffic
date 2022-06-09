"""
Ollivier-Ricci curvatures computed locally using Laplacian [https://arxiv.org/abs/1712.00875]
"""

using Distributed
using Graphs
using JuMP, GLPK
using Graphs.Parallel

∇(f, e, here) = f[findfirst(==(e.dst), here)] - f[findfirst(==(e.src), here)]

function κ(
    G::Tg where Tg <: AbstractGraph, 
    e::Te where Te <: AbstractEdge,
    D::AbstractArray,
    Δ::AbstractArray; 
    optimizer = GLPK.Optimizer
    )

    here = unique( [neighbors(G, e.src)..., neighbors(G, e.dst)...] )
    n = length(here)

    d = view(D, here, here)
    δ = view(Δ, here, here)
    
    model = Model(optimizer)
    @variable(model, f[1:n])

    # f ∈ Lip(1)
    @constraint(model, .- d .≤ [f[x] - f[y] for x ∈ 1:n, y ∈ 1:n] .≤ d) 

    # ∇ₑf = 1
    @constraint(model, ∇(f, e, here) == 1)
    
    # inf ∇ₑΔf
    @objective(model, Min, ∇(δ * f, e, here))

    optimize!(model)

    return objective_value(model)

end

function κ(
    G::AbstractGraph,
    E::Union{Graphs.AbstractEdgeIter, Vector{T} where T <: AbstractEdge} = edges(G);
    optimizer = GLPK.Optimizer,
    parallel = false
    )

    Δ = laplacian_matrix(G)./degree(G)
    @time D = (parallel ? Graphs.Parallel : Graphs).floyd_warshall_shortest_paths(squash(G)).dists

    @time curvatures = (parallel ? pmap : map)(e -> κ(G, e, D, Δ; optimizer = optimizer), E)
    return curvatures

end