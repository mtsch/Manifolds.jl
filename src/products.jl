# TODO: vsak manifold{1} je curve
#       posledično je circle tud
#       mora se convertat št. outputov pred produktom.
#       mogoče rezanje nul na koncu/codimension?

# ScaledManifold
# ReparameterizedManifold

using Random
using LinearAlgebra
using StaticArrays
using ForwardDiff: Dual, partials, value

# = Abstract & helpers = #
struct CurveProduct{D, T<:NTuple{D, Curve}} <: AbstractManifold{D}
    spaces ::T
end
Base.show(io::IO, cp::CurveProduct{D}) where D =
    print(io, join(string.(cp.spaces), " × "))

CurveProduct(spaces::T) where {T<:NTuple{D, Curve}} where D =
    CurveProduct{D, T}(spaces)

function offsetframe(ps::CurveProduct{D}, args::NTuple{D, T}) where {D, T}
    spaces = ps.spaces

    offset = zeros(T, D+3)
    frame  = Matrix{T}(I, (D+3, D+3))

    for i in 1:D
        space = spaces[i]
        scalef = scale(space)
        val, mat = offsetframe(space, args[i])
        offset += frame * vcat(val, zeros(T, D)) *
            (i ≠ 1 ? scalef(args[i-1]) : scalef(zero(T)))
        frame *= idpad(mat, D+3)
    end
    offset, frame
end

(ps::CurveProduct{D})(args::Vararg{T, D}) where {D, T} = first(offsetframe(ps, args))

constant1(_) = 1

#const ScaledCurve = Tuple{Curve, Any}

# Make sure `Curve`s get transformed to `Tuple{Curve, typeof(constant1)}`
#LinearAlgebra.cross(c::Curve, a::Any) = (c, constant1) × a
#LinearAlgebra.cross(a::Any, c::Curve) = a × (c, constant1)
#LinearAlgebra.cross(c1::Curve, c2::Curve) = (c1, constant1) × (c2, constant1)

# All variants of constructing `CurveProduct`s.
LinearAlgebra.cross(c1::Curve, c2::Curve) =
    CurveProduct((c1, c2))
LinearAlgebra.cross(ps::CurveProduct{D}, c::Curve) where D =
    CurveProduct((ps.spaces..., c))
LinearAlgebra.cross(c::Curve, ps::CurveProduct{D}) where D =
    CurveProduct((c, ps.spaces))
LinearAlgebra.cross(p1::CurveProduct, p2::CurveProduct{D}) where D =
    CurveProduct((p1.spaces..., p2.spaces...))

struct PointSpace <: AbstractManifold{0} end
Base.show(io::IO, ::PointSpace) = print(io, "{⋆}")

(ps::PointSpace)() = 0.0
LinearAlgebra.cross(::PointSpace, ::PointSpace) = PointSpace()
LinearAlgebra.cross(m::AbstractManifold, ::PointSpace) = m
LinearAlgebra.cross(::PointSpace, m::AbstractManifold) = m
