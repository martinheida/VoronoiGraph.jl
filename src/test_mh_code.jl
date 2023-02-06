using Revise
using VoronoiGraph
using Revise


function test()
    xs = VoronoiGraph.vecvec(rand(3,100))
    m,S = VoronoiGraph.sysvoronoi(xs,VoronoiGraph.cuboid(3))
    println(length(m.All_Verteces[1])," ",m.All_Verteces[1])
    VoronoiGraph.systematic_refine!(m,VoronoiGraph.vecvec(0.25*rand(3,100)),VoronoiGraph.cuboid(3))

    I = VoronoiGraph.Integrator(m, type=VoronoiGraph.VI_POLYGON, integrand=x->[sum(x),1.0] )
    # auch: VI_GEOMETRY, VI_MONTECARLO
    VoronoiGraph.integrate(I,domain=VoronoiGraph.cuboid(3))

    println(    sum(I.Integral.volumes)    )#length(m.All_Verteces[1])," ",m.All_Verteces[1])
    println(    sum(I.Integral.bulk_integral)    )#length(m.All_Verteces[1])," ",m.All_Verteces[1])
end

test()