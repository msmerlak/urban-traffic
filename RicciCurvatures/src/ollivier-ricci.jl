using Graphs
using JuMP, GLPK #Gurobi


import JuMP.fix
function fix(C::Vector{VariableRef}, V::Vector{T}) where T <: Number
    @assert length(C) == length(V)
    for i in 1:length(C)
        fix(C[i], V[i])
    end
end

function W₁(μ₁::AbstractVector{T}, μ₂::AbstractVector{T}, model::Model) where T <: Real

    fix(model.obj_dict[:marginal1], μ₁)
    fix(model.obj_dict[:marginal1], μ₂)

    optimize!(model)

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

function ollivier_ricci(g::Tg where Tg <: AbstractGraph, edges = Graphs.edges(g); parallel = true, optimizer = GLPK.Optimizer)

    n, _ = size(g)
    distances = floyd_warshall_shortest_paths(g).dists

    model = Model(optimizer)
    
    #set_optimizer_attribute(model, "OutputFlag", 0)
        
    @variable(model, coupling[1:n, 1:n] >= 0)
    @variable(model, marginal1[1:n])
    @variable(model, marginal2[1:n])

    @objective(model, Min, sum(coupling .* distances))

    @constraint(model, normalization, sum(coupling) == 1)

    return [κ(g, e, distances, model) for e ∈ edges]
end




