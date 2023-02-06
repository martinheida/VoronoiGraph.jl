struct RaycastDomainIncircleSkip{T,TT}
    tree::T
    lmesh::Int64
    visited::Vector{Int64}
    ts::Vector{Float64}
    recursive::Bool
    positions::BitVector
    variance_tol::Float64
    break_tol::Float64
    node_tol::Float64
    b_nodes_tol::Float64
    correcting::Bool
    vectors::Matrix{Float64}
    symmetric::Matrix{Float64}
    rhs::Vector{Float64}
    rhs_cg::Vector{Float64}
    ddd::Vector{Float64}
    domain::TT
    rare_events::Vector{Int64}
    dimension::Int64
end

# SRI = search rare index
const SRI_vertex_tolerance_breach = 1
const SRI_vertex_suboptimal_correction = 2
const SRI_vertex_irreparable = 3
const SRI_fraud_boundary_vertex = 4
const SRI_deactivate_boundary = 5
const SRI_activate_mirror = 6
const SRI_irregular_node = 7

const SRI_max = 10


function RaycastDomainIncircleSkip(xs,recursive,variance_tol,break_tol,node_tol,b_tol,correcting,dom,brut)
    lxs=length(xs)
    dim=length(xs[1])
    z1d_1=zeros(Float64,lxs)
    z1d_2=zeros(Float64,dim)
    z1d_3=zeros(Float64,dim)
    z1d_4=zeros(Float64,dim+1)
    z2d_1=zeros(Float64,dim,dim)
    z2d_2=zeros(Float64,dim,dim)
    tree=brut ? BruteTree(xs) : MyTree(xs,length(dom))
    return RaycastDomainIncircleSkip{typeof(tree),typeof(dom)}( tree, lxs, zeros(Int64,lxs+length(dom)+3), z1d_1, recursive, BitVector(zeros(Int8,length(xs))), variance_tol, break_tol, node_tol, b_tol,
                                    correcting, z2d_1, z2d_2, z1d_2, z1d_3, z1d_4, dom, zeros(Int64,SRI_max),dim)
end


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Stepest Decent to find a vertex from thin air

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

function new_vertex!(new_nodes,searcher,r,u,t,nodes,_Cell=nodes[1])
    r2=r+t*u
    if !(r2 in searcher.domain) #&& sig2[end]<=length(xs)
        #println("Hallo1")
        index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in nodes)) 
        if t0<t
            activate_mirror(searcher,_Cell,index)
            for i in 1:(length(new_nodes))
                if !(new_nodes[i] in nodes) 
                    new_nodes[i]=index+searcher.tree.size
                    break
                end
            end
            sort!(new_nodes)
            r2= r+t0*u
        end
    end
    return r2        
end

function adjust_vertex_t(searcher,r,u,t,nodes,_Cell)
        index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in nodes)) 
        if t0<Inf
            activate_mirror(searcher,_Cell,index)
            return sort!(push!(nodes,index+searcher.tree.size)) , t0
        end
    return nodes, t        
end


""" starting at given points, run the ray shooting descent to find vertices """
function descent(xs::Points, searcher::RaycastDomainIncircleSkip, start = 1) 
    sig = [start]
    r = xs[start]
    d = length(r)

        for k in d:-1:1  # find an additional generator for each dimension
            u = randray(xs[sig])
            (tau, t) = raycast(sig, r, u, xs, searcher)
            b=false
            if t==Inf 
                (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start) 
                b=true
            end
            if t == Inf
                u = -u
                (tau, t) = raycast(sig, r, u, xs, searcher)
            end
            if t==Inf 
                (tau,t)=adjust_vertex_t(searcher,r,u,t,sig,start)
                b=true 
            end
            if t == Inf
                error("Could not find a vertex in both directions of current point." *
                    "Consider increasing search range (tmax)")
            end
            if b
                r=r+t*u
                sig=tau
                continue
            end
            r=new_vertex!(tau,searcher,r,u,t,sig,start)
            sig=tau
        end
    sort!(sig)
    r = project(r,searcher.domain)
    r,_ = walkray_correct_vertex(r, sig, searcher, true)
    return (sig, r)
end

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## WALKRAY

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(sig::Sigma, r::Point, xs::Points, searcher::RaycastDomainIncircleSkip, i)
    #global walk_count=walk_count+1
    old = sig[i]
    dim = length(r)
    sig_del = deleteat(sig, i)
    u = randray(xs[sig_del])
    success = true
    if (u' * (xs[sig[i]] - xs[sig_del[1]])) > 0
        u = -u
    end
    repeat = true
    while repeat
        repeat = false
        sig2, t = raycast(sig_del, r, u, xs, searcher, old)
        if t==0.0 
            global SECONDCORRECTIONS+=1
            return sig2, r, u
        end
        #if t<0  println("t=$t ;  ( $sig , $r ) <-> ( $sig2 , $(r + t*u))") end
        if t < Inf
            r2,_ = walkray_correct_vertex(r + t*u, sig2, searcher, false) 
            active = false
            index = 0
            if !(r2 in searcher.domain) #&& sig2[end]<=length(xs)
                index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
                if t0<t
                    active = activate_mirror(searcher,sig_del[1],index)
                    sig2=sort!([sig_del;index+searcher.tree.size])
                    r2 =  r + t0*u 
                end
            end
            r2,success =  walkray_correct_vertex(r2, sig2, searcher, true)
            if active && !success
                deactivate_mirror(searcher,sig_del[1],index)
                return sig2, r2, u, success
            elseif active && success # make sure all nodes of vertex are found on irregular grids, hence need to repeat from beginning with new info. 
                                     # The performance loss on regular grids in 5 dimensions is approximately 0.001
                repeat = true
                continue
            else
                return sig2, r2, u, success
            end
        else
            index, t0 = intersect(searcher.domain,r,u, x-> !((x+searcher.tree.size) in sig_del)) 
            if t0<Inf
                active = activate_mirror(searcher,sig_del[1],index)
                # SRI_fraud_boundary_vertex
                sig2=sort!([sig_del;index+searcher.tree.size])
                r2,success =  walkray_correct_vertex(r + t0*u, sig2, searcher, true) 
                if active && !success
                    deactivate_mirror(searcher,sig_del[1],index)
                elseif active && success # make sure all nodes of vertex are found on irregular grids, hence need to repeat from beginning with new info.
                                     # The performance loss on regular grids in 5 dimensions is approximately 0.001
                    repeat = true
                    continue
                end
                return sig2, r2 , u, success
            else
                return sig_del, r, u, false  # if the vertex has an unbounded ray, return the same vertex
            end
        end
    end
end



function walkray_correct_vertex(_r, sig, searcher, correct_bulk)
    dim = searcher.dimension
    r=_r
    return adjust_boundary_vertex(r, searcher.domain, sig, searcher.lmesh, dim+1, searcher.b_nodes_tol), true
end
   
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## Implementations of two different raycasting search algorithms

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################



######## default raycast

"""
    SysRaycast(xs;recursive=true,variance_tol=1.0E-20,break_tol=1.0E-5,correcting=true,domain=FullSpace(),bruteforce=false)

Initializes the standard searcher for computation of Voronoi meshes.

# Arguments
- `xs::Points`: An array of points, preferably of SVector-type to speed up the algorithm
- `recursive` : When set to 'true' this will cause an iteration on quasi-periodic grids. As long as in some step 'n' verteces are found 
                that primarily belong to vertex 'm<n' restart after each full iteration. This is in praxis only relevant "outside the point cloud"
- `variance-tol`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol` : when the afore mentioned variance is even larger than that (happens only on quasi-periodic grids) this is sign that something goes 
                really wrong. Therefore, the vertex is skipped. Typically happens "outside" the quasi-periodic domain
- `b_nodes_tol`: When a vertex evidently should lie on the boundary but is slightly appart with distance less than `b_nodes_tol`, it will be corrected.
                Otherwise it will be dumped.
- `domain` : When a vertex is found that lies outside of "domain" the algorithm will not look for further neighbor verteces. However, 
            the vertex itself will be stored.
- `bruteforce`: when set to "true" the algorithm will use a BruteTree instead of a KDTree
"""
SysRaycast(xs;recursive=true,variance_tol=1.0E-20,break_tol=1.0E-5,nodes_tol=1.0E-5,b_nodes_tol=1.0E-7,correcting=true,domain=Boundary(),bruteforce=false) = RaycastDomainIncircleSkip(xs,recursive,variance_tol,break_tol,nodes_tol,b_nodes_tol,correcting,domain,bruteforce)



########################################################################################################################################

## Second raycast-method :   RAYCAST Data Structure

########################################################################################################################################

function activate_cell(searcher,_Cell,neigh)
    lxs=searcher.tree.size
    for n in neigh
        if n>lxs 
            activate_mirror(searcher,_Cell,n-lxs) 
        end
    end
end

function activate_mirror(searcher,i,plane)
    if searcher.tree.active[plane] 
        return false 
    end
    searcher.tree.active[plane]=true
    searcher.tree.extended_xs[searcher.tree.size+plane]=reflect(searcher.tree.extended_xs[i],searcher.domain,plane)
    searcher.rare_events[SRI_activate_mirror] += 1
    return true
#    println("node $i : activate $plane <-> $(searcher.tree.size+plane)  ;  $(searcher.tree.extended_xs[i]) <-> $(searcher.tree.extended_xs[searcher.tree.size+plane])")
end

function deactivate_mirror(searcher,i,plane)
    searcher.tree.active[plane] = false
    searcher.rare_events[SRI_deactivate_boundary] += 1
end

struct MyTree{T,TT}
    tree::T
    extended_xs::Vector{TT}
    active::BitVector
    size::Int64
    mirrors::Int64
    function MyTree{T,TT}(t,m,a,s,m2) where {T,TT}
        return new(t,m,a,s,m2)        
    end
    function MyTree(xs,l=0)
        t=KDTree(xs)
        m=append!(copy(xs),Vector{typeof(xs[1])}(undef,l))
        a=BitVector(zeros(Int8,l))
        return MyTree{typeof(t),typeof(xs[1])}(t,m,a,length(xs),l)
    end
end

function _nn(tree::MyTree,x;skip=x->false)
    idx,dists=knn(tree.tree,x,1,false,skip)
    b=length(idx)>0
    index = b ? idx[1] : 0
    dist = b ? dists[1] : Inf64
    lm=tree.mirrors
    if lm==0 return index, dist end
    for i in 1:lm
        ( !tree.active[i] || skip(tree.size+i) ) && continue
        d=norm(x-tree.extended_xs[i+tree.size])
        if d<dist
            index=i+tree.size
            dist=d
        end
    end    
    return index, dist
end



########################################################################################################################################

## raycast-method

########################################################################################################################################


global RAYCAST_ERROR=0::Int64

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastDomainIncircleSkip, old = 0)
    x0 = xs[sig[1]]
    searcher.visited.*=0
    searcher.ts.*=0
    visited = view(searcher.visited,4:length(searcher.visited))
    status = view(searcher.visited,1:3)
    c = maximum(dot(xs[g], u) for g in sig)

    # only consider points on the right side of the hyperplane
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig

    local i, t
        i, t = _nn(searcher.tree, r + u * (u' * (x0-r)), skip=skip)
    t == Inf && return [0], Inf

    # sucessively reduce incircles unless nothing new is found
    k=1
    visited[1]=i
    j=0
    while true
        if (i>length(xs)) || i<=0
            global RAYCAST_ERROR+=1
            return [0], Inf
        end
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        searcher.ts[k]=t
        j, _ = _nn(searcher.tree, r+t*u)
        if j in sig 
            break
        elseif _visited(j,visited)
            searcher.rare_events[9]+=1
            i=j
            break
        else
            i = j
        end
        k+=1
        visited[k]=j
    end
    k=1
    # lv=length(visited)
    while (visited[k]!=i)
        visited[k]=0
        searcher.ts[k]=0
        k+=1
    end
    kk = k
    while visited[kk]!=0
        kk += 1
    end
    if (kk!=k+1)
        searcher.rare_events[SRI_irregular_node] += 1
        tau = sort!(append!(copy(sig),view(visited,k:(kk-1))))
    
        # i == old && (t=0.0) # evtl nicht richtig by degenerate nodes
    
        return tau, t
    else
        visited[k]=0
        searcher.ts[k]=0.0
        # in the end, the above implies that visited[k]!=0 if and only if visited[k] appeared in a cycle and visited[k]!=i 
        tau = sort!([i; sig])
    
        i == old && (t=0.0)
    
        return tau, t
    end
end

function _visited(j,visited)
    k=1
    ret=false
    while (visited[k]!=0)
        if visited[k]==j 
            ret=true 
            break
        end
        k+=1
    end
    return ret
end

