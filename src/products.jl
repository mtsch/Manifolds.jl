# Unit for * and ×.
struct UnitSpace <: Manifold{0, 0} end

Base.one(::Type{<:Manifold}) = UnitSpace()
Base.one(::Manifold) = UnitSpace()

"""
TODO
"""
struct ProductSpace{D, C, N} <: Manifold{D, C}
    spaces::NTuple{N, Manifold}
end

(ps::ProductSpace{D})(args::Vararg{Any, D}) where D =
    ps(promote(args...)...)

function (ps::ProductSpace{D, C, N})(args::Vararg{T, D}) where {T, D, C, N}
    spaces = ps.spaces

    argoffset = 0
    trans = Frame(eye(T, 1), [zero(T)])
    t = zero(T)
    for i in 1:N-1
        d = dim(spaces[i])
        trans *= tnbframe(spaces[i], args[argoffset + (1:d)]..., t)

        t = args[argoffset + (1:d)]
        argoffset += d
    end

    (trans * ps.spaces[end](args[argoffset+1:end]..., scale=t))
end

@inline function tnbframe_undefined_throw(m::Manifold{D}) where D
    T = typeof(m)
    if !method_exists(tnbframe, (T, fill(Float64, D)...))
        error("`tnbframe` not implemented for `$(T)`! ",
              "Use `×` or `cross` instead.")
    end
end

function Base.:*(m1::Manifold{D1, C1},
                 m2::Manifold{D2, C2}) where {D1, D2, C1, C2}
    tnbframe_undefined_throw(m1)
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, 2}((m1, m2))
end

function Base.:*(ps::ProductSpace{D1, C1, N},
                 m::Manifold{D2, C2}) where {D1, D2, C1, C2, N}
    tnbframe_undefined_throw(ps.spaces[end])
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((ps.spaces..., m))
end

function Base.:*(m::Manifold{D1, C1},
                 ps::ProductSpace{D2, C2, N}) where {D1, D2, C1, C2, N}
    tnbframe_undefined_throw(m)
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N+1}((m, ps.spaces...))
end

function Base.:*(ps1::ProductSpace{D1, C1, N1},
                 ps2::ProductSpace{D2, C2, N2}) where {D1, D2, C1, C2, N1, N2}
    tnbframe_undefined_throw(ps1.spaces[end])
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    ProductSpace{D, C, N1+N2}((ps1.spaces..., ps2.spaces...))
end

Base.:*(ps::ProductSpace, us::UnitSpace) = ps
Base.:*(us::UnitSpace, ps::ProductSpace) = ps
Base.:*(m::Manifold, us::UnitSpace) = m
Base.:*(us::UnitSpace, m::Manifold) = m
Base.:*(us1::UnitSpace, us2::UnitSpace) = us1

"""
TODO
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
Base.cross(m::Manifold, us::UnitSpace) = m
Base.cross(us::UnitSpace, m::Manifold) = m
Base.cross(us1::UnitSpace, us2::UnitSpace) = us1
