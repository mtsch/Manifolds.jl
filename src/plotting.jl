# Curve
@recipe function f(curve::Manifold{1, C}, from = 0, to = 1;
                   npoints=1000, projection = eye(3)) where C
    # projection matrix
    X = idpad(projection, (1+C, 3))

    pts = [X'curve(t) for t in linspace(from, to, npoints)]

    linewidth --> 2
    color     --> :orange
    legend    --> false

    xs = getindex.(pts, 1)
    ys = getindex.(pts, 2)
    zs = getindex.(pts, 3)

    xs, ys, zs
end

# Surface
@recipe function f(surface::Manifold{2, C}, from = [0,0], to = [1,1];
                   npoints = 100, projection = eye(3)) where C
    # projection matrix
    X = idpad(projection, (2+C, 3))

    xs = zeros(npoints, npoints)
    ys = zeros(npoints, npoints)
    zs = zeros(npoints, npoints)
    for (i, θ) in enumerate(linspace(from[1], to[1], npoints)),
        (j, φ) in enumerate(linspace(from[2], to[2], npoints))
        xs[i, j], ys[i, j], zs[i, j] = X'surface(θ, φ)
    end
    seriestype := :surface
    xs, ys, zs
end

# TODO: unsatisfying
@recipe function f(s::Manifold{3}, from = [0,0,0], to = [1,1,1];
                   npoints = 10)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]

    for (i, θ) in enumerate(linspace(from[1], to[1], npoints)),
        (j, φ) in enumerate(linspace(from[2], to[2], npoints)),
        (k, ψ) in enumerate(linspace(from[3], to[3], npoints))
        x, y, z = s(θ, φ, ψ)[1:3]
        push!(xs, x)
        push!(ys, y)
        push!(zs, z)
    end
    seriestype := :scatter
    markersize --> 0.5
    xs, ys, zs
end

# Transformation
@recipe function f(fr::Frame)
    p = idpad(fr.translation, 3)
    B = 0.5idpad(fr.basis_change, 3) .+ p

    for (i, col, coord) in zip(1:3, [:red, :green, :blue], [:x, :y, :z])
        @series begin
            color := col
            label := coord
            linewidth --> 2
            [p[1], B[1, i]], [p[2], B[2, i]], [p[3], B[3, i]]
        end
    end
end

# 3d points.
@userplot Points3d
@recipe function f(p::Points3d; projection = eye(3))
    if length(p.args) != 1 || !(typeof(first(p.args)) <: AbstractMatrix)
        error("points3d is expecting a single matrix argument. " *
              "Got: $(typeof(p.args))")
    end

    @assert size(projection, 2) == 3
    pts = first(p.args)
    dim = size(pts, 1)

    # projection matrix
    X = idpad(projection, (dim, 3))

    pts = X'pts

    seriestype := :scatter
    markersize --> 0.5
    pts[1, :], pts[2, :], pts[3, :]
end
