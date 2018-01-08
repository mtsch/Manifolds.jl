"""
    Manifold{D, C}

A manifold of dimension `D` and codimension `C`. A `Manifold` must be callable
with `D` arguments and should return a vector of length `D + C`. The call should
also include a `scale` keyword argument that is used to scale its overall size.
"""
abstract type Manifold{Dimension, Codimension} end

dim(::Manifold{D}) where D = D
codim(::Manifold{D, C}) where {D, C} = C
ambientdim(::Manifold{D, C}) where {D, C} = D + C

Base.copy(man::Manifold) = man

Base.rand(man::Manifold, n; scale = nothing, noise = 0) =
    rand(Base.GLOBAL_RNG, man, n, scale = scale, noise = noise)

function Base.rand(rng::AbstractRNG, man::Manifold{D, C}, n = 1;
                   scale = nothing, noise = 0) where {D, C}
    res = Matrix{Float64}(D + C, n)
    for i in 1:n
        hiss = noise * randn(D + C)
        if scale == nothing
            res[:, i] = man(rand(rng, D)...) .+ hiss
        else
            res[:, i] = man(rand(rng, D)..., scale=scale) .+ hiss
        end
    end
    res
end

"""
    Interval(f::Union{Funtcion, Number})

An interval scaled by `f`.
"""
struct Interval <: Manifold{1, 0}
    scale::Function
end

Interval(l::Number = 1.0) = Interval(_ -> l)

(int::Interval)(t; scale = 0.0) = SVector(int.scale(scale) * t)

"""
    NSphere{N}(r::Union{Funtcion, Number})

A N-Sphere with radius `r`. Circles and 3-Spheres can be constructed with
`Circle(r)` and `Sphere(r)`, respectively.
"""
struct NSphere{N} <: Manifold{N, 1}
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
struct Knot{T} <: Manifold{1, 2}
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
