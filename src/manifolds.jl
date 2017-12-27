# ============================================================================ #
abstract type Manifold{Dimension, Codimension} end

dim(::Manifold{D}) where D = D
codim(::Manifold{D, C}) where {D, C} = C
ambientdim(::Manifold{D, C}) where {D, C} = D + C

# Allow arithmetic.
struct UnitSpace <: Manifold{0, 0} end

Base.oneunit(::Type{<:Manifold}) = UnitSpace()
Base.one(::Type{<:Manifold}) = UnitSpace()
Base.oneunit(::Manifold) = UnitSpace()
Base.one(::Manifold) = UnitSpace()
Base.copy(man::Manifold) = man

# TODO: specialize for Manifold{1, 2} -- binormal always [0,0,1]?
# TODO: lifting might be wrong with high dimensional spaces.
# TODO: TransformedManifold - A*m where A is a matrix?
"""
    tnbframe(man::Manifold, ts...; scale=nothing)

Return a point and coordinate frame perpendicular to the manifold `man`
evaluated at `man(ts..., scale=scale)` wrapped in a `RigidTransformation`.
"""
function tnbframe(man::Manifold{1, C}, t::T, scale=nothing) where {C, T}
    C > 2 && error("Cannot calculate the tnbframe!") # TODO: figure this out - it might work.

    t_dual = Dual{:d1}(Dual{:d2}(t, one(T)), one(T))
    if scale == nothing
        val = idpad(man(t_dual), 3)
    else
        val = idpad(man(t_dual, scale=scale), 3)
    end

    tangent  = normalize(value.(partials.(val, 1)))
    normal   = partials.(partials.(val, 1), 1)
    binormal = normalize(tangent × normal)
    normal   = binormal × tangent

    RigidTransformation(SMatrix{3,3}([tangent normal binormal]),
                        value.(value.(val)))
end

function tnbframe(man::Manifold{1, 0}, t::T, scale=nothing) where T
    if scale == nothing
        val = man(t)[1]
    else
        val = man(t, scale = scale)[1]
    end
    RigidTransformation(eye(SMatrix{3, 3, T}), @SVector([val[1], 0, 0]))
end

function tnbframe(man::Manifold{2}, t1, t2, scale=nothing)
    if scale == nothing
        val = idpad(man(t1, t2), 3)
        J   = jacobian(v -> man(v...), [t1, t2])
    else
        val = idpad(man(t1, t2, scale = scale), 3)
        J   = jacobian(v -> man(v[1:end-1]..., scale=v[end]), [t1, t2, scale])
    end

    tangent  = normalize(J[:, 1])
    normal   = J[:, 2]
    binormal = normalize(tangent × normal)
    normal   = binormal × tangent

    RigidTransformation(SMatrix{3, 3}([tangent normal binormal]), val)
end

# TODO: support other number types -- add type parameter to Manifold*.
Base.rand(man::Manifold, n; scale=nothing, noise=0) =
    rand(Base.GLOBAL_RNG, man, n, scale=scale, noise=noise)

function Base.rand(rng::AbstractRNG, man::Manifold{D, C}, n=1;
                   scale=nothing, noise=0) where {D, C}
    res = Matrix{Float64}(D+C, n)
    for i in 1:n
        hiss = noise * randn(D+C)
        if scale == nothing
            res[:, i] = man(rand(rng, D)...) .+ hiss
        else
            res[:, i] = man(rand(rng, D)..., scale=scale) .+ hiss
        end
    end
    res
end

"""
TODO.
"""
struct Interval <: Manifold{1, 0}
    scale::Function
end

Interval(l::Real=1.0) = Interval(_ -> l)

(int::Interval)(t; scale=0.0) = SVector(int.scale(scale) * t)

"""
TODO.
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

function Base.rand(rng::AbstractRNG, s::NSphere{N}, n=1;
                   scale=nothing, noise=0) where N

    res = Matrix{Float64}(N+1, n)
    for i in 1:n
        hiss = noise * randn(N+1)
        res[:, i] = s.scale(scale) .* normalize!(randn(rng, N+1)) .+ hiss
    end
    res
end

"""
TODO.
"""
# TODO: come up with nice names for parameters.
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
        trans = lift(trans(tnbframe(spaces[i], args[argoffset + (1:d)]..., t)), d)

        t = args[argoffset + (1:d)]
        argoffset += d
    end

    trans(ps.spaces[end](args[argoffset+1:end]..., scale=t))[1:D+C]
end

@inline function tnbframe_exists_throw(m::Manifold{D}) where D
    T = typeof(m)
    if !method_exists(tnbframe, (T, fill(Float64, D)...))
        error("`tnbframe` not implemented for `$(T)`! ",
              "Use `×` or `cross` instead.")
    end
end

function Base.:*(m1::Manifold{D1, C1},
                 m2::Manifold{D2, C2}) where {D1, D2, C1, C2}
    tnbframe_exists_throw(m1)
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, 2}((m1, m2))
end

function Base.:*(ps::ProductSpace{D1, C1, N},
                 m::Manifold{D2, C2}) where {D1, D2, C1, C2, N}
    tnbframe_exists_throw(ps.spaces[end])
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((ps.spaces..., m))
end

function Base.:*(m::Manifold{D1, C1},
                 ps::ProductSpace{D2, C2, N}) where {D1, D2, C1, C2, N}
    tnbframe_exists_throw(m)
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((m, ps.spaces...))
end

function Base.:*(ps1::ProductSpace{D1, C1, N1},
                 ps2::ProductSpace{D2, C2, N2}) where {D1, D2, C1, C2, N1, N2}
    tnbframe_exists_throw(ps1.spaces[end])
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N1+N2}((ps1.spaces..., ps2.spaces...))
end

Base.:*(ps::ProductSpace, us::UnitSpace) = ps
Base.:*(us::UnitSpace, ps::ProductSpace) = ps

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
    for i in 1:N
        d = dim(spaces[i])
        c = codim(spaces[i])

        currargs = args[argoffset + (1:d)]
        res[resoffset + (1:d+c)] .= spaces[i](currargs...)

        argoffset += d
        resoffset += d+c
    end

    SVector{D+C}(res)
end

function Base.rand(rng::AbstractRNG,
                   cs::CartesianSpace{D, C, N}, n=1) where {D, C, N}
    spaces = cs.spaces
    res = Matrix{Float64}(D+C, n)

    resoffset = 0
    for i in 1:N
        a = ambientdim(spaces[i])

        res[resoffset + (1:a), :] .= rand(rng, spaces[i], n)
        resoffset += a
    end
    res
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
    CartesianSpace{D1+D2, C1+C2, N1+N2}((cs1.spaces..., cs2.spaces...))
end

Base.cross(cs::CartesianSpace, us::UnitSpace) = cs
Base.cross(us::UnitSpace, cs::CartesianSpace) = cs
