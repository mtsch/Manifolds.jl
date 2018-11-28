# Abstract, helpers ====================================================================== #
"""
    AbstractManifold{D, C}

A manifold of dimension `D` and codimension `C`.

Should be callable with `D` arguments in [0, 1] and should return a `SVector{D+C}`.
"""
abstract type AbstractManifold{D, C} end

dim(::AbstractManifold{D}) where D = D
codim(::AbstractManifold{<:Any, C}) where C = C

Base.rand(man::AbstractManifold) = rand(Random.GLOBAL_RNG, man)
Base.rand(rng::AbstractRNG, man::AbstractManifold{D}) where D =
    man(rand(Random.GLOBAL_RNG, D)...)
Base.rand(man::AbstractManifold, n::Integer) = rand(Random.GLOBAL_RNG, man, n)
Base.rand(rng::AbstractRNG, man::AbstractManifold{D}, n::Integer) where D =
    [rand(rng, man) for _ in 1:n]

# Point ================================================================================== #
"""
    PointSpace()

Single point, unit for `×`.
"""
struct PointSpace <: AbstractManifold{0, 1} end
Base.show(io::IO, ::PointSpace) = print(io, "{⋆}")

(ps::PointSpace)() = SVector(0.0)
LinearAlgebra.cross(::PointSpace, ::PointSpace) = PointSpace()
LinearAlgebra.cross(m::AbstractManifold, ::PointSpace) = m
LinearAlgebra.cross(::PointSpace, m::AbstractManifold) = m

# Spheres ================================================================================ #
"""
    Sphere{D}

A D-Sphere.
"""
struct Sphere{D} <: AbstractManifold{D, 1}
end

Sphere(D) = Sphere{D}()
Base.show(io::IO, ::Sphere{D}) where D = print(io, "𝕊^$D")

function (::Sphere{D})(varargs::Vararg{T, D}) where {D, T}
    args = 2 .* varargs
    SVector((prod(sinpi.(args[1:i-1])) * cospi(args[i]) for i in 1:D)...,
            prod(sinpi.(args[1:end])))
end

Base.rand(rng::AbstractRNG, ::Sphere{D}) where D =
    normalize(SVector{D + 1, Float64}(randn(rng, D + 1)))

# Ball =================================================================================== #
"""
    Ball{D}

A D-Ball.
"""
struct Ball{D} <: AbstractManifold{D, 0}
end

Ball(D) = Ball{D}()
Base.show(io::IO, ::Ball{D}) where D = print(io, "𝔹^$D")

function (::Ball{D})(varargs::Vararg{T, D}) where {D, T}
    args = 2 .* varargs
    (args[1]/2)^(1/D) .*
        SVector((prod(sinpi.(args[2:i-1])) * cospi(args[i]) for i in 2:D)...,
                prod(sinpi.(args[2:end])))
end

function Base.rand(rng::AbstractRNG, ::Ball{D}) where D
    res = @SVector fill(2.0, D)
    while norm(res) > 1
        res = SVector{D, Float64}(2rand(rng, D) .- 1)
    end
    res
end

# Cube =================================================================================== #
"""
    Cube{D}

A D-Cube.
"""
struct Cube{D} <: AbstractManifold{D, 0}
end

Cube(D) = Cube{D}()
Base.show(io::IO, ::Cube{D}) where D = print(io, "𝕀", D ≠ 1 ? "^$D" : "")
(::Cube{D})(args::Vararg{T, D}) where {D, T} = SVector(args)

# ParametricManifold ===================================================================== #
struct ParametricManifold{D, C, F} <: AbstractManifold{D, C}
    fun::F
end

ParametricManifold{D, C}(f::F) where {D, C, F} =
    ParametricManifold{D, C, F}(f)
Base.show(io::IO, pm::ParametricManifold{D, C}) where {D, C} =
    print(io, "ParametricManifold{$D, $C}($(pm.fun))")
(pm::ParametricManifold{D, C})(args::Vararg{T, D}) where {D, C, T} =
    SVector{D+C, T}(pm.fun(args...))
