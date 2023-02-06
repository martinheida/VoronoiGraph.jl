module VoronoiGraph

using LinearAlgebra
using NearestNeighbors
using Polyhedra
using ProgressMeter
using SparseArrays
using StaticArrays
using RecipesBase

include("voronoi.jl")
include("raycast.jl")
include("volume.jl")
include("plot.jl")
include("montecarlo.jl")

include("mesh.jl")
include("boundary.jl")
include("integral.jl")
include("sysvoronoi.jl")
include("meshrefine.jl")
include("sysraycast.jl")
include("integrate.jl")
include("mcintegrator.jl")
include("Leibnitzrule.jl")
include("polyintegrator.jl")

export voronoi, voronoi_random
export volumes
export mc_volumes, mc_integrate

for i in 1:10
    precompile(voronoi, (Vector{SVector{i, Float64}},))
end

end
