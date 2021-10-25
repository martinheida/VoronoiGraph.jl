using SpecialFunctions


"""
    mc_volume(i, xs, nmc, searcher)

Estimate the area and volume of the `i`-th Voronoi cell from the Voronoi Diagram generated by `xs`
by a Monte Carlo estimate from random rays
"""
function mc_volume(i, xs, nmc=1000, searcher = SearchIncircle(100, KDTree(xs)))

    x = xs[i]
    d = length(x)

    V = 0.
    A = spzeros(length(xs))

    for iter in 1:nmc
        u = normalize(randn(d))
        (j, t) = raycast([i], x, u, xs, searcher)

        V += t^d

        t == Inf && continue
        j = pop!(setdiff(j, i))
        normal = normalize(xs[j] - xs[i])
        A[j] += t ^ (d-1) / abs(dot(normal, u))
    end

    V *= π^(d/2) / gamma(d/2 + 1) / nmc
    A *= 2 * π^(d/2) / gamma(d/2) / nmc

    A, V
end

"""
    mc_volumes(xs::Points, nmc=1000)

Estimate the areas and volumes of the Voronoi Cells generated by `xs` using `nmc` Monte Carlo samples.

# Returns
- `SparseMatrix`: the areas of the common boundaries of two cells
- `Vector`: the volumes of each cell
"""
function mc_volumes(xs::Points, nmc=1000)
    V = zeros(length(xs))
    A = spzeros(length(xs), length(xs))

    searcher = SearchIncircle(100, KDTree(xs))

    for i in 1:length(xs)
        a, v = mc_volume(i, xs, nmc, searcher)
        V[i] = v
        A[:, i] = a
    end
    #A = (A + A') / 2
    A, V
end
