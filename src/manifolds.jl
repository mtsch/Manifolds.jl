# ============================================================================ #
abstract type Manifold{Dimension, Codimension} end
# TODO: change E to codimension for spheres.
# TODO: make spheres, knots etc into types.

dim(::Manifold{D}) where D = D
codim(::Manifold{D, C}) where {D, C} = C
ambientdim(::Manifold{D, C}) where {D, C} = D + C

# TODO: specialize for Manifold{1, 2}?
# TODO: build lift into frame
function frame(man::Manifold{1, C}, t::T, tα=nothing) where {C, T}
    C > 2 && error("Cannot calculate the frame!") # TODO: figure this out.

    τ = Dual{:d1}(Dual{:d2}(t, one(T)), one(T))
    if tα == nothing
        val = idpad(man(τ), 3)
    else
        val = idpad(man(τ, scale=tα), 3)
    end

    n = normalize(value.(partials.(val, 1)))
    b = partials.(partials.(val, 1), 1)
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), value.(value.(val)))
end

function frame(pc::Manifold{1, 0}, t::T, tα=nothing) where T
    if tα == nothing
        RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t)[1],0,0]))
    else
        RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t, scale=tα)[1],0,0]))
    end
end

function frame(ps::Manifold{2}, θ, φ, tα=nothing)
    if tα == nothing
        pt = idpad(ps(θ, φ), 3)
        J  = jacobian(v -> ps(v...), [θ, φ])
    else
        pt = idpad(ps(θ, φ, scale = tα), 3)
        J  = jacobian(v -> ps(v[1:end-1]..., scale=v[end]), [θ, φ, tα])
    end

    n = normalize(J[:, 1])
    b = J[:, 2]
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), pt)
end

"""
"""
struct NSphere{N} <: Manifold{N, 1}
    scale::Function
end

NSphere{N}(r::Real = 1.0) where N = NSphere{N}(_ -> r)

function (s::NSphere{N})(varargs::Vararg{T, N}; scale=0) where {N, T}
    res  = Vector{T}(N+1)
    args = 2π .* varargs
    r    = s.scale(scale)
    r .* SVector((prod(sin.(args[1:i-1])) * cos(args[i]) for i in 1:N)...,
                 prod(sin.(args[1:end])))
end

"""
"""
struct Interval <: Manifold{1, 0}
    scale::Function
end

Interval(l::Real=1.0) = Interval(_ -> l)

(int::Interval)(t; scale=0.0) = SVector(int.scale(scale) * t)

"""
"""
struct Knot{T} <: Manifold{1, 2}
    p::T
    q::T
    n::T
    m::T
    h::T
    scale::Function
end

function Knot(;p=2.0, q=3.0, n=1.0, m=1.5, h=1.0, scale=1.0)
    p, q, n, m, h = promote(p,q,n,m,h)
    T = typeof(p)
    if isa(scale, Function)
        Knot{T}(p, q, n, m, h, scale)
    else
        Knot{T}(p, q, n, m, h, _ -> scale)
    end
end

function (kn::Knot)(t; scale=0)
    p = kn.p
    q = kn.q
    n = kn.n
    m = kn.m
    h = kn.h
    kn.scale(scale) .* @SVector([m*cos(2π*p*t) + n*cos((2-q)*2π*t),
                                 m*sin(2π*p*t) + n*sin((2-q)*2π*t),
                                 h * sin(q*2π*t)])
end

# TODO: cartesian products (ala clifford torus)
"""
"""
struct ProductSpace{D, C, N} <: Manifold{D, C}
    spaces::NTuple{N, Manifold}
end

(ps::ProductSpace{D})(args::Vararg{Any, D}) where D =
    ps(promote(args...)...)
function (ps::ProductSpace{D, C, N})(args::Vararg{T, D}) where {T, D, C, N}
    spaces = ps.spaces

    argoffset = 0
    trans = id(1, T)
    t = zero(T)
    for i in 1:N-1
        d = dim(spaces[i])
        trans = lift(trans(frame(spaces[i], args[argoffset + (1:d)]..., t)), d)

        t = args[argoffset + (1:d)]
        argoffset += d
    end

    trans(ps.spaces[end](args[argoffset+1:end]..., scale=t))[1:D+C]
end

function Base.:*(m1::Manifold{D1, C1},
                 m2::Manifold{D2, C2}) where {D1, D2, C1, C2}
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, 2}((m1, m2))
end

function Base.:*(ps::ProductSpace{D1, C1, N},
                 m::Manifold{D2, C2}) where {D1, D2, C1, C2, N}
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((ps.spaces..., m))
end

function Base.:*(m::Manifold{D1, C1},
                 ps::ProductSpace{D2, C2, N}) where {D1, D2, C1, C2, N}
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((m, ps.spaces...))
end

function Base.:*(ps1::ProductSpace{D1, C1, N1},
                 ps2::ProductSpace{D2, C2, N2}) where {D1, D2, C1, C2, N1, N2}
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N1+N2}((ps1.spaces..., ps1.spaces...))
end

"""
"""
struct CartesianSpace{D, C, N} <: Manifold{D, C}
    spaces::NTuple{N, Manifold}
end

(cs::CartesianSpace{D})(args::Vararg{Any, D}) where D =
    cs(promote(args...)...)
function (cs::CartesianSpace{D, C, N})(args::Vararg{T, D}) where {T, D, C, N}
    spaces = cs.spaces
    res = Vector{T}(D+C)

    argoffset = 0
    resoffset = 0
    lastargs  = 0:0
    for i in 1:N
        d = dim(spaces[i])
        c = codim(spaces[i])
        currargs = args[argoffset + (1:d)]

        res[resoffset + (1:d+c)] .= spaces[i](currargs..., scale=lastargs)

        argoffset += d
        resoffset += d+c
        lastargs = currargs
    end

    SVector{D+C}(res)
end

function Base.cross(m1::Manifold{D1, C1},
                    m2::Manifold{D2, C2}) where {D1,D2,C1,C2}
    CartesianSpace{D1+D2, C1+C2, 2}((m1, m2))
end

function Base.cross(cs::CartesianSpace{D1, C1, N},
                    m::Manifold{D2, C2}) where {D1,D2,C1,C2,N}
    CartesianSpace{D1+D2, C1+C2, N+1}((cs.spaces..., m))
end

function Base.cross(m::Manifold{D1, C1},
                    cs::CartesianSpace{D2, C2, N}) where {D1,D2,C1,C2,N}
    CartesianSpace{D1+D2, C1+C2, N+1}((m, cs.spaces...))
end

function Base.cross(cs1::CartesianSpace{D1, C1, N1},
                    cs2::CartesianSpace{D2, C2, N2}) where {D1,D2,C1,C2,N1,N2}
    CartesianSpace{D1+D2, C1+C2, N1+N2}((cs1.spaces..., cs1.spaces...))
end

#=
# ============================================================================ #
interval(α=1.0) =
    ParametricCurve(1, t -> SVector(t), α)

circle(r=1.0) =
    ParametricCurve(2, t -> @SVector([cos(2π*t), sin(2π*t)]), r)

knot(p=2, q=3, n=1, m=1.5, h=1, α=1) =
    ParametricCurve(3, t -> @SVector([m*cos(2π*p*t) + n*cos((2-q)*2π*t),
                                      m*sin(2π*p*t) + n*sin((2-q)*2π*t),
                                      h * sin(q*2π*t)]), α)

torus(R, r, α=1.0) =
    ParametricSurface(3, (θ, φ) -> @SVector([(R + r*cos(2π*θ))*cos(2π*φ),
                                             (R + r*cos(2π*θ))*sin(2π*φ),
                                             r*sin(2π*θ)]), α)

sphere(α=1.0) =
    ParametricSurface(3, (θ, φ) -> @SVector([cos(2π*φ)*sin(π*θ),
                                             sin(2π*φ)*sin(π*θ),
                                             cos(π*θ)]), α)

kleinbottle(R, P, ϵ=0, α=1.0) =
    ParametricSurface(4, (θ, φ) ->
                      @SVector([R*(cos(π*θ)*cos(2π*φ) - sin(π*θ)*sin(4π*φ)),
                                R*(sin(π*θ)*cos(2π*φ) - cos(π*θ)*sin(4π*φ)),
                                P*cos(2π*θ) * (1 + ϵ*sin(2π*φ)),
                                P*sin(2π*θ) * (1 + ϵ*sin(2π*φ))]), α)

# ============================================================================ #
struct ParametricCurve{E} <: Manifold{1, E}
    fun::Function
    scale::Function
end

# Constant scale
ParametricCurve(edim::Int, f::Function, scale = 1.0) =
    ParametricCurve{edim}(f, _ -> scale)

# Function scale
ParametricCurve(edim::Int, f::Function, scalefun::Function) =
    ParametricCurve{edim}(f, scalefun)

function (pc::ParametricCurve)(t, s = 0)
    t1, s1 = promote(t, s)
    pc.fun(t1) .* pc.scale(s1)
end

# ============================================================================ #
struct ParametricSurface{E} <: Manifold{2, E}
    fun::Function
    scale::Function
end

# Constant scale
ParametricSurface(edim::Int, f::Function, scale = 1.0) =
    ParametricSurface{edim}(f, _ -> scale)

# Function scale
ParametricSurface(edim::Int, f::Function, scalefun::Function) =
    ParametricSurface{edim}(f, scalefun)

function (ps::ParametricSurface)(θ, φ, s = 0)
    t1, t2, s1 = promote(θ, φ, s)
    ps.fun(t1, t2) .* ps.scale(s1)
end

=#
