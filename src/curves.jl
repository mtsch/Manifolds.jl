struct ParametricCurve{F<:Function, T<:Real}
    fun::F
    trans::RigidTransform
    dim::Int
    scale::T
end

ParametricCurve(f::F, d::Int=1, α::T=1) where {F<:Function, T<:Real} =
    ParametricCurve{F, T}(f, id(d, T), d, α)

(pc::ParametricCurve)(t) =
    pc.trans(pc.fun(t) .* pc.scale)

(t::RigidTransform)(pc::ParametricCurve) =
    ParametricCurve(pc.fun, t ∘ pc.trans, pc.dim, pc.scale)

function frenetframe(pc::ParametricCurve, t::T) where T #TODO NE T
    τ = hyper(t, 1.0, 1.0, 0.0)
    pt = pc.fun(τ) * pc.scale
    pt3 = idpad(pt, 3)

    n = normalize(eps1.(pt3))
    b = eps1eps2.(pt3)
    t = normalize(n × b)
    b = -(n × t)

    pc.trans(SMatrix{3, 3, T}([n b t]), false), pc.trans(real.(pt))
end

function lift(t::RigidTransform{D}, n=1) where D
    if n < 1
        t
    else
        t ∘ change_basis([zeros(n, D) I;
                          I zeros(D, n)])
    end
end

function lift(pc::ParametricCurve, n=1)
    if n < 1
        pc
    else
        ParametricCurve(pc.fun, lift(pc.trans, n), pc.dim, pc.scale)
    end
end

circlefun(t) = @SVector [cos(2π*t), sin(2π*t)]
trefoilfun(t) = @SVector [sin(2π*t) + 2sin(4π*t),
                          cos(2π*t) - 2cos(4π*t),
                         -sin(6π*t)]
circle(r) = ParametricCurve(circlefun, 2, r)
trefoil(r) = ParametricCurve(trefoilfun, 3, r)

struct ParametricManifold{N}
    curves::NTuple{N, ParametricCurve}
end

function (pm::ParametricManifold{N})(args::Vararg{T, N}) where {T, N}
    curves = pm.curves

    fm, pt = frenetframe(curves[1], args[1])
    trans = lift(RigidTransform(fm, pt))
    for i in 2:N
        fm, pt = frenetframe(trans(curves[i]), args[i])
        trans = lift(RigidTransform(fm, pt), i)
    end

    # Cut off excess dimensions.
    pt[1:(curves[end].dim + N-1)]
end

function Base.cross(args::Vararg{ParametricCurve, N}) where N
    ParametricManifold{N}(args)
end

# ======== GARBAGE ========= #
if false

    trf = (translate([1,1,1]) ∘ rotate_x(π/2))(trefoil(1.0))
    plt = plot(trf)
    for t in linspace(0, 1, 100)
        #point = idpad(trf(t), 3)
        frame, point = frenetframe(trf, t)
        mframe = frame .+ point
        colors = [:red, :green, :blue]
        # Frenet frame.
#       for i in 1:3
#           plot!(plt,
#                 [point[1], mframe[1, i]],
#                 [point[2], mframe[2, i]],
#                 [point[3], mframe[3, i]],
#                 color=colors[i])
#       end
        c = lift(RigidTransform(frame, point))(circle(0.1))
        plot!(c)
    end
    plt


    @time begin
        torus = cross(trefoil(3.0), circle(1.0))
        mat = zeros(3, 10000)
        for i in 1:10000
            mat[:, i] .= torus(rand(), rand())
        end
    end
    scatter(mat[1,:], mat[2,:], mat[3,:], markersize=0.5)

    # TODO to recipe
    torus = cross(trefoil(3.0), circle(1.0))
    xs = zeros(100, 100)
    ys = zeros(100, 100)
    zs = zeros(100, 100)
    for (i, θ) in enumerate(linspace(0,1,100)),
        (j, φ) in enumerate(linspace(0,1,100))
        x,y,z = torus(θ, φ)
        xs[i, j] = x
        ys[i, j] = y
        zs[i, j] = z
    end
    surface(xs', ys', zs')

end
