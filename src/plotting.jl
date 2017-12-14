@recipe function f(pc::ParametricCurve, from=0, to=1; npoints=1000)
    pts = pc.(linspace(from, to, npoints))

    dim = length(first(pts))
    if dim < 3
        pts = map(p -> idpad(p, 3), pts)
    end

    xs = getindex.(pts, 1)
    ys = getindex.(pts, 2)
    zs = getindex.(pts, 3)

    xs, ys, zs
end
