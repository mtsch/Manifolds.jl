# remove codim
"""
    AbstractAbstractManifold{D, C}

A manifold of dimension `D` and codimension `C`. An `AbstractAbstractManifold` must be
callable with `D` arguments and should return a tuple of length `D + C`.
"""
abstract type AbstractAbstractManifold{D, C} end

dim(::AbstractManifold{D}) where D = D
codim(::AbstractManifold{D, C}) where {D, C} = C
ambientdim(::AbstractManifold{D, C}) where {D, C} = D + C

Base.copy(man::AbstractManifold) = man

Base.rand(man::AbstractManifold, n) =
    rand(Base.GLOBAL_RNG, man, n)

function Base.rand(rng::AbstractRNG, man::AbstractManifold{D, A}, n = 1) where {D, A}
    [SVector{D + C, T}(man(rand(rng, D)...)) for _ in 1:n]
end

struct ParametricManifold{D, A, F} <: AbstractManifold{D, A}
    f::F
end

function curve(f, codim)

end

"""
    NSphere{N}(r::Union{Funtcion, Number})

A N-Sphere with radius `r`. Circles and 3-Spheres can be constructed with
`Circle(r)` and `Sphere(r)`, respectively.
"""
struct NSphere{N} <: AbstractManifold{N, 1}
    r::Float64
end

NSphere{N}(r = 1.0) where N = NSphere{N}(Float64(r))

function (s::NSphere{N})(varargs::Vararg{T, N}) where {N, T}
    args = 2π .* varargs
    s.r .* SVector((prod(sin.(args[1:i-1])) * cos(args[i]) for i in 1:N)...,
                   prod(sin.(args[1:end])))
end

function Base.rand(rng::AbstractRNG, s::NSphere{N}, n = 1) where N

    res = Vector{SVector{N+1, Float64}}(undef, n)
    for i in 1:n
        res[i] = s.r .* normalize!(@SVector randn(rng, N + 1))
    end
    res
end


#=
"""
    Interval(f::Union{Funtcion, Number})

An interval scaled by `f`.
"""
struct Interval <: AbstractManifold{1, 0}
    scale::Function
end

Interval(l::Number = 1.0) = Interval(_ -> l)

(int::Interval)(t; scale = 0.0) = SVector(int.scale(scale) * t)

"""
    NSphere{N}(r::Union{Funtcion, Number})

A N-Sphere with radius `r`. Circles and 3-Spheres can be constructed with
`Circle(r)` and `Sphere(r)`, respectively.
"""
struct NSphere{N} <: AbstractManifold{N, 1}
    scale::Function
end

NSphere{N}(r::Number = 1.0) where N = NSphere{N}(_ -> r)

function (s::NSphere{N})(varargs::Vararg{T, N}; scale = 0) where {N, T}
    args = 2π .* varargs
    rad  = s.scale(scale)
    rad .* SVector((prod(sin.(args[1:i-1])) * cos(args[i]) for i in 1:N)...,
                   prod(sin.(args[1:end])))
end

function Base.rand(rng::AbstractRNG, s::NSphere{N}, n = 1;
                   scale = nothing, noise = 0) where N

    res = Matrix{Float64}(N + 1, n)
    for i in 1:n
        hiss = noise * randn(N + 1)
        res[:, i] = s.scale(scale) .* normalize!(randn(rng, N + 1)) .+ hiss
    end
    res
end

"""
    Knot{T}(; p = 2, q = 3, n = 1, m = 1.5, h = 1, scale = 1)

A knot parameterized by `p`, `q`, `n`, `m` and `h` and scaled by `scale`:

    x = m*cos(2π*p*t) + n*cos((2-q)*2π*t)
    y = m*sin(2π*p*t) + n*sin((2-q)*2π*t)
    z = h * sin(q*2π*t)])
"""
# TODO: come up with nice names for parameters.
struct Knot{T} <: AbstractManifold{1, 2}
    p::T
    q::T
    n::T
    m::T
    h::T
    scale::Function
end

function Knot(;p = 2.0, q = 3.0, n = 1.0, m = 1.5, h = 1.0, scale = 1.0)
    p, q, n, m, h = promote(p, q, n, m, h)
    T = typeof(p)
    if scale isa Function
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

=#
