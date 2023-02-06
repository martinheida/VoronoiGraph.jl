####################################################################################################################################

## Most simple integrator to calculate the verteces of the mesh and potentially the neighbors

####################################################################################################################################


struct Geometry_Integrator
    Integral::Voronoi_Integral
    function Geometry_Integrator(mesh::Voronoi_MESH,neigh=false)
        N=(Vector{Int64})[]
        if neigh
            l=length(mesh)
            emptyint=Int64[]
            N=Vector{Vector{Int64}}(undef,l)
            for i in 1:l N[i]=copy(emptyint) end
        end
        return new(Voronoi_Integral{typeof(mesh.nodes[1])}(N,[],[],[],[],mesh))
    end
    function Geometry_Integrator(points::Points,neigh=false)
        return Geometry_Integrator(Voronoi_MESH(points),neigh)
    end
    function Geometry_Integrator(Inte::Voronoi_Integral)
        return new(Inte)
    end
end

function copy(I::Geometry_Integrator)
    return Geometry_Integrator(copy(I.Integral))    
end

function integrate(xs,c,a,b,s,I::Geometry_Integrator)
end

function integrate(Integrator::Geometry_Integrator; domain=nothing, relevant=1:2, modified=1:2) 
    for i in 1:(length(Integrator.Integral.neighbors))
        Integrator.Integral.neighbors[i]=neighbors_of_cell(i,Integrator.Integral.MESH.All_Verteces[i],Integrator.Integral.MESH.Buffer_Verteces[i])
    end
end

function prototype_bulk(Integrator::Geometry_Integrator)
    return Float64[]
end

function prototype_interface(Integrator::Geometry_Integrator)
    return Float64[]
end

####################################################################################################################################

## Two different functions to delete entries in arrays

####################################################################################################################################




function my_deleteat(sig,i,dim)
    newsig=Array{Int64}(undef,dim) #zeros(Int64,dim)
    k=1
    while k<i
        newsig[k]=sig[k] #push!(newsig,sig[k])
        k+=1
    end
    while k<=dim
        newsig[k]=sig[k+1] #push!(newsig,sig[k+1]) #
        k+=1
    end
    return newsig
end



#############################################################################################

# Here follows everything special about systematic Voronoi search

#############################################################################################

function sysvoronoi(xs::Points,boundary=Boundary(); searcher=SysRaycast(xs,domain=boundary),Iter=1:length(xs)) 
    return sysvoronoi(Voronoi_MESH(xs),boundary,searcher=searcher,Iter=Iter)
end
 
function sysvoronoi(mesh::Voronoi_MESH,boundary=Boundary(); Iter=1:(length(mesh.nodes)), searcher=SysRaycast(mesh.nodes,domain=boundary)) 

    #mesh=Integrator.Integral.MESH
    xs=searcher.tree.extended_xs
    dimension=length(xs[1])

    boundary = mesh.boundary_Verteces
    
    edgecount_global = Dict{Vector{Int64},Char}()

    l=length(xs)
    if l==0 || l<=dimension
        error("There are not enough points to create a Voronoi tessellation")
    end

    TODO = collect(Iter) 
    repeat=true
    iteration_count=1
    TODO_count=length(TODO)
    new_verteces=0::Int64
    while repeat && iteration_count<=4
        repeat=false
        k=1
        while k<=length(TODO) # iterate:
            i=TODO[k]
            TODO[k]=0
            k+=1
            i==0 && continue
            # Expolaration of neighbor verteces, given the verteces of all "previous" cells, i.e. all smaller relevant values of i
            new_verteces+=systematic_explore_cell(xs,i,mesh,edgecount_global,boundary,searcher)
        end
        if searcher.recursive
            count=1
            for j in 1:searcher.tree.size
                if searcher.positions[j] && (j in Iter) 
                    TODO[count]=j
                    count+=1
                    repeat=true
                end
            end
            TODO_count=count-1
            searcher.positions.*=0
        end 
        iteration_count+=1
        #break
    end
    println(searcher.rare_events)
    return mesh, searcher
end

##############################################################################################################################

## Core functions of the geometry part

##############################################################################################################################

global NO_VERTEX=0::Int64
global VER_VAR=0.0::Float64
global CORRECTIONS=0::Int64
global SUCCESSFUL=0::Int64
global FIRSTCORRECTIONS=0::Int64
global SECONDCORRECTIONS=0::Int64

function systematic_explore_cell(xs::Points,_Cell,mesh,edgecount_local,boundary,searcher,offset=0)
    new_verteces=0::Int64
    dimension=length(xs[1])
    lmesh = length(mesh)
    #load all known vertices of the current cell
    verts = mesh.Buffer_Verteces[_Cell]
    allverts = mesh.All_Verteces[_Cell]
    searcher.tree.active.*=0
    neigh=neighbors_of_cell(_Cell,allverts,verts)
    activate_cell(searcher,_Cell,neigh)
    #initialize container to collect edges of _Cell
    empty!(edgecount_local)
    # create an empty list where newly found verteces will be stored intermittently
    queue = EmptyDictOfType(Int64[]=>xs[1]) #copy(allverts)
    lxs=length(xs)
    if length(verts)!=0
        for (sig,r) in verts
            if sig[2]==_Cell # && sig[dimension+1]<=lxs
                edge=deleteat(sig,1)
                edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
            end
        end
        # modify here if all local edges shall be collected!!!
    end

    if length(allverts)!=0
        for (sig,r) in allverts
            start = 2 #sig[dimension+1]<=lxs ? 2 : dimension+1 
            for j in start:(dimension+1) # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
                edge = my_deleteat(sig, j, dimension)
                edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
            end
        end
    end

    if length(verts)==0 &&  length(allverts)==0 #i.e. if length(queue)==0
        sig=[0]
        r=xs[1]
        k=0
        while (sig[1]<_Cell) # in quasi-periodic media, sig[1]<_Cell with positive probability
            sig, r = descent(xs,searcher,_Cell)
            if !(r in searcher.domain)
                println("so ein mist: $sig -> $r")
            end
            push!(queue, sig=>r)
            for j in 1:length(sig)
                edge = deleteat(sig, j)
                c=get(edgecount_local, edge, '0') + 1
                if c>='3'
                    k+=1
                    break
                end
                edgecount_local[edge] = c
            end
            k>=10 && break
        end
        #if sig[1]!=_Cell println("problem at $_Cell: found $sig") end
    end


    for (sig,r) in allverts
        #print("  $sig ")
        
        if (sig[2]<_Cell)# || sig[end]>lxs) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
        #if !(r in searcher.domain) continue end
        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            println("at $_Cell: should never happen - explore_cell 1: $sig")
            systematic_explore_vertex(xs,sig,r,_Cell,1,edgecount_local,verts,queue,allverts,boundary,searcher) 
                    # call method once for a replacement of the first entry, all other verteces with sig[1] and
                    # sig[2]=_Cell  are already known by former iterations 
            continue
        end
        # otherwise, we know sig[1]=_Cell. All other nodes of "sig" have not been visited or replaced yet. 
        # So we fix this entry and start our search the other "d" directions
        l=length(sig)-offset
        
        for i in 2:l
            systematic_explore_vertex(xs,sig,r,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
        end
    end

    for (sig,r) in verts
        #print("  $sig ")
        
        if (sig[2]<_Cell)# || sig[end]>lxs) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
        #if !(r in searcher.domain) continue end
        if (sig[1]!=_Cell)# && !get(edgecount_global, deleteat(sig, i), '0') != '2') # in case the vertex was found in an earlier step, but the second entry is already _Cell
            
            systematic_explore_vertex(xs,sig,r,_Cell,1,edgecount_local,verts,queue,allverts,boundary,searcher) 
                    # call method once for a replacement of the first entry, all other verteces with sig[1] and
                    # sig[2]=_Cell  are already known by former iterations 
            continue
        end
    end

    while length(queue) > 0
        new_verteces+=1
        (sig,rr) = pop!(queue)
        r=rr
        #println("  Queue: $sig")
        sig[end]>lxs && continue # this seems to be a relic which probably never happens
        if (sig[2]<_Cell) # in case the vertex was already part in at least two iterations, everything about this vertex is known.
            continue
        end
 
        if (sig[1]!=_Cell)
            haskey(all_verts[sig[1]],sig) && continue # if we do not use this line, we run into trouble in quasi-periodic meshes
                    # in high dimensions. This should not be, maybe it is a strange Julia-Problem. However, the additional safty-call
                    # will not cost much, as this case is rather rare... If you do not believe me, erase this line and activate the following
            #println("here $sig, $_Cell, $(haskey(all_verts[sig[1]],sig)), $(haskey(buffer_verts[sig[2]],sig))")
            push!(mesh,sig=>r,lmesh)
            searcher.positions[sig[1]]=true
            searcher.positions[sig[2]]=true
            continue
        end
        # otherwise, we know sig[1]=_Cell. All other nodes of "sig" have not been visited or replaced yet. 
        # So we fix this entry and start our search the other "d" directions
        #push!(allverts,sig=>r)
        push!(mesh,sig=>r,lmesh)
        l=length(sig)-offset

        for i in 2:l
            systematic_explore_vertex(xs,sig,r,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
        end

    end

    return new_verteces
end


function systematic_explore_vertex(xs,sig,R,_Cell,i,edgecount_local,verts,queue,allverts,boundary,searcher)
    oldnode=sig[i]
    dimension=length(xs[1])
    #print(" $i($(sig[i])): ")
    edge=my_deleteat(sig, i, dimension) # current edge to walk along
    if get(edgecount_local, edge, '0') >= '2'  #if edge explored, cancel routine 
        return
    end
    new_sig, r,u, success = walkray(sig, R, xs, searcher, i) # provide missing node "j" of new vertex and its coordinate "r" 
                                                    # together with edge orientation 'u'
    if length(sig) > length(new_sig) #if oldnode==newnode then we found a boundary element and we can cancel 
        push!(boundary, new_sig=>boundary_vertex(R,u,sig[i]))
        return
    end

    if !success
        return
    end

    #in all other cases, we have a potentially new relevant vertex
    newvertex=new_sig
    if !haskey(allverts, newvertex) && !haskey(queue, newvertex) #in case we really have new vertex ....
        push!(queue, newvertex => r)  # put it to the queue
        #push!(allverts, newvertex => r)
        for j in 2:length(new_sig) # Allways: new_sig[1]=_Cell and we only need edges with edge[1]=_Cell
            edge = my_deleteat(new_sig, j, dimension)
            edgecount_local[edge] = get(edgecount_local, edge, '0') + 1
        end
    end
    return
end

 