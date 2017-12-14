abstract type Manifold{D, E} end

dim(::Manifold{D}) where D = D
embeddedin(::Manifold{D, E}) where {D, E} = E

struct ParametricCurve{E} <: Manifold{1, E}
    fun::Function
    scale::Function
end

# Constant scale
ParametricCurve(edim::Int, f, scale = 1.0) =
    ParametricCurve{edim}(f, _ -> scale)

# Function scale
ParametricCurve(edim::Int, f, scalefun::Function) =
    ParametricCurve{edim}(f, scalefun)

# TODO: specialize for ParametricCurve{2}
function frame(pc::ParametricCurve, t::T, α=0) where T
    τ = hyper(t, one(T), one(T), zero(T))
    #pt = pc.fun(τ) .* pc.scale
    pt  = pc(τ, α)
    pt3 = idpad(pt, 3)

    n = normalize(eps1.(pt3))
    b = eps1eps2.(pt3)
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), real.(pt3))
end

function frame(pc::ParametricCurve{1}, t::T, α=0) where T
    RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t, α)[1],0,0]))
end

(pc::ParametricCurve)(t, s = 0.0) =
    pc.fun(t) .* pc.scale(s)

# ============================================================================ #
struct ProductSpace{D, E, N} <: Manifold{D, E}
    spaces::NTuple{N, Manifold}
end

function (ps::ProductSpace{D, E, N})(args::Vararg{T, D}) where {T, D, E, N}
    spaces = ps.spaces

    argoffset = 0
    trans = id(1, T)
    α = zero(T)
    for i in 1:N-1
        d = dim(spaces[i])
        trans = lift(trans(frame(spaces[i], args[argoffset + (1:d)]..., α)), d)

        α = args[argoffset + 1]
        argoffset += d
    end

    trans(spaces[end](args[argoffset+1:end]..., α))[1:E]

    #=
    fm, pt = frenetframe(curves[1], args[1])
    trans = lift(RigidTransform(fm, pt))
    for i in 2:N
        fm, pt = frenetframe(trans(curves[i]), args[i])
        trans = lift(RigidTransform(fm, pt), i)
    end

    # Cut off excess dimensions.
    pt[1:(curves[end].dim + N-1)]
    =#
end

function Base.:*(m1::Manifold{D1, E1},
                 m2::Manifold{D2, E2}) where {D1, D2, E1, E2}
    ProductSpace{D1+D2, E1+D2, 2}((m1, m2))
end

function Base.:*(ps::ProductSpace{D1, E1, N},
                 m::Manifold{D2, E2}) where {D1, D2, E1, E2, N}
    ProductSpace{D1+D2, E1+D2, N+1}((ps.spaces..., m))
end

function Base.:*(m::Manifold{D1, E1},
                 ps::ProductSpace{D2, E2, N}) where {D1, D2, E1, E2, N}
    ProductSpace{D1+D2, E1+D2, N+1}((m, ps.spaces...))
end

function Base.:*(ps1::ProductSpace{D1, E1, N1},
                 ps2::ProductSpace{D2, E2, N2}) where {D1, D2, E1, E2, N1, N2}
    ProductSpace{D1+D2, E1+D2, N1+N2}((ps1.spaces..., ps1.spaces...))
end





# == EXAMPLES == #
circle(r=1.0) =
    ParametricCurve(2, t -> @SVector([cos(2π*t), sin(2π*t)]), r)

trefoil(α=1.0) =
    ParametricCurve(3, t -> @SVector([sin(2π*t) + 2sin(4π*t),
                                      cos(2π*t) - 2cos(4π*t),
                                     -sin(6π*t)]), α)

interval(α=1.0) =
    ParametricCurve(1, t -> SVector(t), α)

function torusknot(p::Int, q::Int, r=2, α=1.0)
    gcd(p, q) ≠ 1 && warn("`p` and `q` are not relatively prime!")
    ParametricCurve(3, t -> @SVector([r*cos(2π*p*t),
                                      r*sin(2π*p*t),
                                     -sin(2π*q*t)]), α)
end

# == GARBAGE == #
if false

    torus = trefoil() * circle(.2)
    plt = plot(torus)
    for i in 1:10
        θ, φ = rand(2)#0., 1.
        p = torus(θ, φ)
        J = ForwardDiff.jacobian(t -> torus(t[1], t[2]), [θ, φ])
        scatter!([p[1]], [p[2]], [p[3]], markersize=2)
        j1 = normalize(J[:, 1])
        j2 = normalize(J[:, 2])
        j3 = -(j1 × j2)
        for (ji, c) in zip([j1,j2,j3], [:red, :green, :blue])
            plot!(plt,
                  p[1] + [0,ji[1]],
                  p[2] + [0,ji[2]],
                  p[3] + [0,ji[3]],
                  color = c,
                  linewidth=2)
        end
    end
    plt

end
