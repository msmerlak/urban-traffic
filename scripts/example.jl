using Distributed
addprocs(5)

@everywhere using DrWatson
@everywhere @quickactivate
@everywhere using RicciCurvatures

using Graphs

g = SimpleGraph(4);
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)

κ(g; parallel = true)