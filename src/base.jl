using Random
using LinearAlgebra
using ForwardDiff: Dual, value, partials

# Abstract, helpers ====================================================================== #
abstract type AbstractManifold{D} end

Base.rand(man::AbstractManifold, n::Integer) = rand(Random.GLOBAL_RNG, man, n)
Base.rand(rng::AbstractRNG, man::AbstractManifold{D}, n::Integer) where D =
    [man(rand(rng, D)...) for _ in 1:n]

"""
    offsetframe(man::AbstractManifold{D}, args::NTuple{D, T})

Return a matrix that makes vectors perpendicular to `man(args...)` and the value of
`man(args...)`.
"""
offsetframe

function idpad(M::AbstractMatrix{T}, d) where T
    n, m = size(M)
    @boundscheck d ≥ n && d ≥ m ||
        throw(ArgumentError("`d` must be larger than the size of `M`"))
    res = zeros(T, d, d)
    res[1:n, 1:m] .= M
    for i in 1:min(d-n, d-m)
        res[n+i, m+i] = 1
    end
    res
end

# Curves ================================================================================= #
struct ParametricCurve{F} <: AbstractManifold{1}
    f::F
end
Base.show(io::IO, c::ParametricCurve) = print(io, "ParametricCurve($(c.f))")

function (c::ParametricCurve)(t::T) where T
    r = c.f(t)
    dim = length(r)
    @boundscheck dim ≤ 3 || error("Invalid curve function result: $r")
    SVector{3, T}(r..., ntuple(zero, 3 - dim)...)
end

const Curve = AbstractManifold{1}

function offsetframe(c::Curve, t::T) where T
    val = c(Dual{:d1}(Dual{:d2}(t, one(T)), one(T)))
    if length(val) < 3
        val = append!([val...], fill(zero(eltype(val)), 3 - length(val)))
    elseif length(val) > 3
        throw(ArgumentError("Curve has codimension $(length(val) - 1)"))
    end

    offset  = value.(value.(val))
    tangent = normalize(value.(partials.(val, 1)))
    normal  = partials.(partials.(val, 1), 1)
    if all(iszero, normal)
        offset, T[0 0;
                  1 0;
                  0 1]
    else
        binormal = normalize(tangent × normal)
        normal   = binormal × tangent

        offset, [normal binormal]
    end
end

# Spheres, Balls ========================================================================= #
"""
    Sphere{D}

A D-Sphere.
"""
struct Sphere{D} <: AbstractManifold{D}
end

Sphere(D) = Sphere{D}()
Base.show(io::IO, ::Sphere{D}) where D = print(io, "S^$D")

function (::Sphere{D})(varargs::Vararg{T, D}) where {D, T}
    args = 2 .* varargs
    SVector((prod(sinpi.(args[1:i-1])) * cospi(args[i]) for i in 1:D)...,
            prod(sinpi.(args[1:end])))
end

Base.rand(rng::AbstractRNG, ::Sphere{D}, n::Integer) where D =
    normalize.([SVector{D + 1, Float64}(randn(rng, D + 1)) for _ in 1:n])

"""
    Ball{D}

A D-Ball.
"""
struct Ball{D} <: AbstractManifold{D}
end

Ball(D) = Ball{D}()
Base.show(io::IO, ::Ball{D}) where D = print(io, "B^$D")

function (::Ball{D})(varargs::Vararg{T, D}) where {D, T}
    args = 2 .* varargs[2:end]
    varargs[1] .* SVector((prod(sinpi.(args[1:i-1])) * cospi(args[i]) for i in 1:D-1)...,
                          prod(sinpi.(args[1:end])))
end

function Base.rand(rng::AbstractRNG, ::Ball{D}, n::Integer) where D
    res = SVector{D, Float64}[]
    while n ≠ 0
        new = filter!(p -> norm(p) ≤ 1, [SVector{D, Float64}(2rand(rng, D) .- 1) for _ in 1:n])
        n = n - length(new)
        append!(res, new)
    end
    res
end
