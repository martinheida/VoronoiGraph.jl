const VI_MIN=0
const VI_POLYGON=2
const VI_MONTECARLO=3
const VI_GEOMETRY=4

const VI_MAX=6



function Integrator(mesh::Voronoi_MESH;type=VI_GEOMETRY,integrand=nothing,bulk_integrand=nothing,interface_integrand=nothing,mc_accurate=(1000,100,20),integral=nothing)
    bi = typeof(bulk_integrand)!=Nothing
    ii = typeof(interface_integrand)!=Nothing
    f = nothing
    fb = nothing
    fi = nothing
    if bi && ii 
        f = x->vcat(bulk_integrand(x),interface_integrand(x))
        fb = bulk_integrand
        fi = interface_integrand 
    elseif bi 
        f = x->bulk_integrand(x)
        fb = f 
    elseif ii
        f = x->interface_integrand(x)
        fi = f
    else
        f = integrand
        fb = f 
        fi = f
    end
    if type==VI_POLYGON
        return Polygon_Integrator(mesh,f,true) 
#    elseif type==VI_HEURISTIC
#        return Heuristic_Integrator(mesh,f,true) 
    elseif type==VI_MONTECARLO
        i,b,r=mc_accurate
        return Montecarlo_Integrator(mesh, b=fb, i=fi, nmc_bulk=b, nmc_interface=i, recycle=r)
    elseif type==VI_GEOMETRY
        return Geometry_Integrator(mesh,true) # let the integrator also calculate the neighbors of the cell
    else
        error("$type is not a valid Integrator. We only allow values between $VI_TEST and $VI_MAX")
    end
end

###############################################################################################################

## stores geometric data for integration

###############################################################################################################

struct IntegrateData
    extended_xs::Points
    domain::Boundary
    size::Int64
    active::BitVector
    float_vec_buffer::Vector{Float64}
    float_vec_vec_buffer::Vector{Vector{Float64}}
    dimension::Int64
    function IntegrateData(xs,dom)
        l=length(dom)
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return new(m,dom,length(xs),a,Float64[],(Vector{Float64})[],length(xs[1]))
    end
end

function activate_data_cell(tree::IntegrateData,_Cell,neigh)
    tree.active .*= 0
    lxs=tree.size
    for n in neigh
        if n>lxs 
            plane=n-lxs
            tree.active[plane]=true
            tree.extended_xs[lxs+plane]=reflect(tree.extended_xs[_Cell],tree.domain,plane)
        end
    end
end


###############################################################################################################

## actual integration method

###############################################################################################################
"""
For each implemented Integrator type this method shall be overwritten. 
In particular, the passage of calculate and iterate might be modified according to the needs of the respective class.
See also Polygon_Integrator and Montecarlo_Integrator for reference.
"""
function integrate(Integrator; domain=Boundary(), relevant=1:(length(Integrator.Integral)+length(domain)), modified=1:(length(Integrator.Integral))) 
    _integrate(Integrator; domain=domain, calculate=relevant, iterate=modified) 
end

"""
Iterates integrate_cell over all elements of iterate. 
It thereby passes the information on whether volume, areas, bulk- or surface integrals shall be calculated.
"""
function _integrate(Integrator; domain=Boundary(), calculate=1:(length(Integrator.Integral)+length(domain)), iterate=1:(length(Integrator.Integral)), 
                    I_data=nothing) 
    TODO=collect(iterate)
    Integral=Integrator.Integral
    data = typeof(I_data)!=Nothing ? I_data : IntegrateData(Integral.MESH.nodes,domain)

    vol=length(Integral.volumes)>0
    ar=length(Integral.area)>0
    bulk=length(Integral.bulk_integral)>0
    inter=length(Integral.interface_integral)>0
    TODO_count=length(TODO)
    max_string_i = length(string(iterate[end], base=10)) 
    max_string_todo = length(string(TODO_count, base=10)) 

    for k in 1:TODO_count # initialize and array of length "length(xs)" to locally store verteces of cells
        integrate_cell(vol,ar,bulk,inter,TODO[k],iterate, calculate, data,Integrator)
    end
    return Integrator
end

"""
adjusts the entries of the Integrator.Integral variable: 
It sorts the entries according to the modified order of neighbors and fills up gaps and deletes entries for neighbors that are gone. 
afterwards it calls the true integration function that is provided by the Integrator.
"""
function integrate_cell(vol::Bool,ar::Bool,bulk::Bool,inter::Bool,  _Cell, iterate, calculate, data, Integrator)
    new_neighbors=neighbors_of_cell(_Cell,Integrator.Integral.MESH.All_Verteces[_Cell],Integrator.Integral.MESH.Buffer_Verteces[_Cell])
    old_neighbors=Integrator.Integral.neighbors[_Cell]
    #println(old_neighbors)
    #println(new_neighbors)
    I=Integrator.Integral
    #isdefined(I.area,_Cell) && println("$(I.area[_Cell])")
    proto_bulk=prototype_bulk(Integrator)
    proto_interface=prototype_interface(Integrator)
    if (length(old_neighbors)>0)
        #print(" ho  ")
        if bulk && (!(isdefined(I.bulk_integral,_Cell)) || length(I.bulk_integral[_Cell])!=length(proto_bulk))
            I.bulk_integral[_Cell]=proto_bulk
        end
        if inter && !(isdefined(I.interface_integral,_Cell))
            I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors))
            for i in 1:(length(old_neighbors)) 
                (I.interface_integral[_Cell])[i]=copy(proto_interface) 
            end
        end
        knn = 0
        for n in new_neighbors
            knn += (n in old_neighbors) ? 0 : 1
        end
        if (knn>0) 
            a_neighbors = zeros(Int64,knn)
            a_areas = zeros(Int64,knn)
            n_interface = Vector{Vector{Float64}}(undef,inter ? knn : 0)
            for i in 1:(inter ? knn : 0)
                n_interface[i]=copy(proto_interface) 
            end
            knn2 = 1
            for n in new_neighbors
                if !(n in old_neighbors)
                    a_neighbors[knn2] = n
                    knn2 += 1
                end
            end
            areas = I.area[_Cell]
            append!(old_neighbors,a_neighbors)
            append!(areas,a_areas)
            inter && append!(I.interface_integral[_Cell],n_interface)
            for k in 1:length(old_neighbors)
                if !(old_neighbors[k] in new_neighbors)
                    old_neighbors[k] = length(data.extended_xs)+data.size
                end
            end
            quicksort!(old_neighbors, ar ? areas : old_neighbors, inter ? I.interface_integral[_Cell] : old_neighbors)
            lnn = length(new_neighbors)
            resize!(old_neighbors,lnn)
            resize!(areas,lnn)
            inter && resize!(I.interface_integral[_Cell],lnn)
        end
    else
        old_neighbors=new_neighbors
        Integrator.Integral.neighbors[_Cell]=new_neighbors
        vol && (I.volumes[_Cell]=0)
        ar && (I.area[_Cell]=zeros(Float64,length(old_neighbors)))
        bulk && (I.bulk_integral[_Cell]=prototype_bulk(Integrator))
        inter && (I.interface_integral[_Cell]=Vector{Vector{Float64}}(undef,length(old_neighbors)))
        inter && (for i in 1:(length(old_neighbors)) 
            (I.interface_integral[_Cell])[i]=copy(proto_interface) 
        end)
    end
    activate_data_cell(data,_Cell,old_neighbors)
    dfvb=data.float_vec_buffer
    dfvvb=data.float_vec_vec_buffer
    #println("$(I.area[_Cell])")
    V=integrate(old_neighbors,_Cell,iterate, calculate, data,Integrator, ar ? I.area[_Cell] : dfvb , bulk ? I.bulk_integral[_Cell] : dfvb , inter ? I.interface_integral[_Cell] : dfvvb)
    if (vol)
        I.volumes[_Cell]=V
    end
end
