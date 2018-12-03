@recipe function f(curve::AbstractManifold{1, C}; npoints=1000) where C
    pts = [curve(t) for t in range(0.0, 0.1, length=npoints)]

    linewidth --> 2
    color     --> :orange
    legend    --> false

    xs = get.(pts, 1, 0.0)
    ys = get.(pts, 2, 0.0)
    zs = get.(pts, 3, 0.0)

    xs, ys, zs
end

@recipe function f(surface::AbstractManifold{2, C}; npoints = 100) where C
    rng = range(0.0, 1.0, length=npoints)

    xs = zeros(npoints, npoints)
    ys = zeros(npoints, npoints)
    zs = zeros(npoints, npoints)
    ws = zeros(npoints, npoints)
    for (i, θ) in enumerate(rng), (j, φ) in enumerate(rng)
        v = surface(θ, φ)
        xs[i, j] = get(v, 1, 0.0)
        ys[i, j] = get(v, 2, 0.0)
        zs[i, j] = get(v, 3, 0.0)
        ws[i, j] = get(v, 4, 0.0)
    end
    seriestype := :surface
    if !iszero(ws)
        surfacecolor := ws
    end

    xs, ys, zs
end

@userplot Points3d
@recipe function f(p::Points3d)
    if length(p.args) != 1 || !(typeof(first(p.args)) <: AbstractVector{<:AbstractVector})
        error("points3d is expecting a vector of vectors. " *
              "Got: $(typeof(p.args))")
    end
    pts = first(p.args)

    seriestype := :scatter
    markersize --> 0.5
    get.(pts, 1, 0.0), get.(pts, 2, 0.0), get.(pts, 3, 0.0)
end
