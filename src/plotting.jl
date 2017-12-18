# Curve
@recipe function f(c::Manifold{1, C}, from=0, to=1; npoints=1000) where C
    pts = c.(linspace(from, to, npoints))

    if C ≤ 2
        pts = map(p -> idpad(p, 3), pts)
    end

    xs = getindex.(pts, 1)
    ys = getindex.(pts, 2)
    zs = getindex.(pts, 3)

    xs, ys, zs
end

# Surface
@recipe function f(s::Manifold{2, C}, from=[0,0], to=[1,1];
                   npoints = 100, project_to = eye(3)) where C
    @assert size(project_to, 2) == 3

    # projection matrix
    Δ = C+2-size(project_to, 1)
    if Δ > 0
        X = [project_to; zeros(Δ, 3)]
    else
        X = project_to
    end

    xs = zeros(npoints, npoints)
    ys = zeros(npoints, npoints)
    zs = zeros(npoints, npoints)
    for (i, θ) in enumerate(linspace(from[1], to[1], npoints)),
        (j, φ) in enumerate(linspace(from[2], to[2], npoints))
        x, y, z = (X'idpad(s(θ, φ), 3))[1:3]
        xs[i, j] = x
        ys[i, j] = y
        zs[i, j] = z
    end
    seriestype := :surface
    xs, ys, zs
end

# TODO
@recipe function f(s::Manifold{3}, from=[0,0,0], to=[1,1,1];
                   npoints = 10) #where N

    #N= 3
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
@recipe function f(rt::RigidTransformation)
    p = point(rt)
    B = 0.5basis(rt) .+ p

    for (i, col, coord) in zip(1:3, [:red, :green, :blue], [:x, :y, :z])
        @series begin
            color := col
            label := coord
            linewidth --> 2
            [p[1], B[1, i]], [p[2], B[2, i]], [p[3], B[3, i]]
        end
    end
end
