
using Graphs
using JuMP, GLPK #Gurobi
using LinearAlgebra:normalize
using Distributed:pmap
using Suppressor

# import JuMP.fix
# function fix(C::Vector{VariableRef}, V::Vector{T}) where T <: Number
#     @assert length(C) == length(V)
#     for i in 1:length(C)
#         fix(C[i], V[i])
#     end
# end

function W₁(μ₁::AbstractVector{T}, μ₂::AbstractVector{T}, model::Model) where T <: Real


    if haskey(model, :marginal1)
        delete(model, model[:marginal1])
        unregister(model, :marginal1)
    end

    if haskey(model, :marginal2)
        delete(model, model[:marginal2])
        unregister(model, :marginal2)
    end

    @constraint(model, marginal1, vec(mapslices(sum, model[:coupling]; dims = 1)) .== μ₁)
    @constraint(model, marginal2, vec(mapslices(sum, model[:coupling]; dims = 2)) .== μ₂)

    @suppress optimize!(model)


    return objective_value(model)
end

function m(g::Tg, x::Int) where Tg <: AbstractGraph
    return Vector( normalize(weights(g)[x, :], 1) ) 
end

function κ(
    g::Tg where Tg <: AbstractGraph, 
    e::Te where Te <: AbstractEdge, 
    d::AbstractMatrix{T} where T <: Real,
    model::Model
    )

    return 1 - W₁(m(g, e.src), m(g, e.dst), model)/d[e.src, e.dst]

end

function ollivier_ricci(g::Tg where Tg <: AbstractGraph, edges = Graphs.edges(g); parallel = false, optimizer = GLPK.Optimizer)

    n, _ = size(g)
    distances = floyd_warshall_shortest_paths(g).dists

    model = Model(optimizer)

    @variable(model, coupling[1:n, 1:n] >= 0)
    @objective(model, Min, sum(coupling .* distances))
    #@constraint(model, normalization, sum(coupling) == 1)

    return (parallel ? pmap : map)(e -> κ(g, e, distances, model), edges)
end
