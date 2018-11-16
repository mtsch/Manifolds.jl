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
struct CurveProduct{D, T<:NTuple{D, Curve}, S<:NTuple{D, Any}} <: AbstractManifold{D}
    spaces ::T
    scales ::S
end
function Base.show(io::IO, cp::CurveProduct{D}) where D
    print(io, "CurveProduct of $D spaces")
end

CurveProduct(spaces::T, scales::S) where {T<:NTuple{D, Curve}, S<:NTuple{D, Any}} where D =
    CurveProduct{D, T, S}(spaces, scales)

function offsetframe(ps::CurveProduct{D}, args::NTuple{D, T}) where {D, T}
    spaces = ps.spaces
    scales = ps.scales

    offset = zeros(T, D+3)
    frame  = Matrix{T}(I, (D+3, D+3))

    for i in 1:D
        space = spaces[i]
        scale = scales[i]
        val, mat = offsetframe(space, args[i])
        offset += frame * vcat(val, zeros(T, D)) *
            (i ≠ 1 ? scale(args[i-1]) : scale(zero(T)))
        frame *= idpad(mat, D+3)
    end
    offset, frame
end

(ps::CurveProduct{D})(args::Vararg{T, D}) where {D, T} = first(offsetframe(ps, args))

constant1(_) = 1

const ScaledCurve = Tuple{Curve, Any}

# Make sure `Curve`s get transformed to `Tuple{Curve, typeof(constant1)}`
LinearAlgebra.cross(c::Curve, a::Any) = (c, constant1) × a
LinearAlgebra.cross(a::Any, c::Curve) = a × (c, constant1)
LinearAlgebra.cross(c1::Curve, c2::Curve) = (c1, constant1) × (c2, constant1)

# All variants of constructing `CurveProduct`s.
LinearAlgebra.cross((c1, scale1)::ScaledCurve, (c2, scale2)::ScaledCurve) =
    CurveProduct((c1, c2), (scale1, scale2))
LinearAlgebra.cross(ps::CurveProduct{D}, (c, scale)::ScaledCurve) where D =
    CurveProduct((ps.spaces..., c), (ps.scales..., scale))
LinearAlgebra.cross((c, scale)::ScaledCurve, ps::CurveProduct{D}) where D =
    CurveProduct((c, ps.spaces), (scale, ps.scales...))
LinearAlgebra.cross(p1::CurveProduct, p2::CurveProduct{D}) where D =
    CurveProduct((p1.spaces..., p2.spaces...), (p1.scales..., p2.scales...))

struct PointSpace <: AbstractManifold{0} end
Base.show(io::IO, ::PointSpace) = print(io, "{⋆}")

(ps::PointSpace)() = 0.0
LinearAlgebra.cross(::PointSpace, ::PointSpace) = PointSpace()
LinearAlgebra.cross(m::AbstractManifold, ::PointSpace) = m
LinearAlgebra.cross(::PointSpace, m::AbstractManifold) = m
