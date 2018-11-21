# Frames ================================================================================= #
"""
    offsetframe(man::AbstractManifold{D, C}, args::NTuple{D, T})

Return a matrix that makes vectors perpendicular to `man(args...)` when multiplied with
`*ₚ` and the value of `man(args...)`.

Should return a `SVector{D+C}` and a matrix.
"""
offsetframe(man::AbstractManifold{D, C}, t::Vararg{T, D}) where {D, C, T} =
    man(t...), vcat(zeros(T, C, 1), one(T))

function offsetframe(c::AbstractManifold{1, C}, t::T) where {C, T}
    # This only works for codimension up to 2.
    if C > 2
        return invoke(offsetframe, Tuple{AbstractManifold{1, C}, T}, c, t)
    end

    val = c(Dual{:d1}(Dual{:d2}(t, one(T)), one(T)))
    if length(val) < 3
        val = append!([val...], fill(zero(eltype(val)), 3 - length(val)))
    elseif length(val) > 3
        # This should never happen, except if C wrong.
        throw(ArgumentError("Curve has codimension $(length(val) - 1)"))
    end

    offset  = value.(value.(val))
    tangent = value.(partials.(val, 1))
    normal  = partials.(partials.(val, 1), 1)
    if iszero(normal × tangent)
        normal = !iszero(tangent × @SVector[0, 0, 1]) ?
            tangent × @SVector[0, 0, 1] : tangent × @SVector[0, 1, 0]
    end
    binormal = normalize(tangent × normal)
    normal   = normalize(binormal × tangent)

    offset[1:1+C], [normal binormal]
end


"""
    u +ₚ v

Promote vectors `u` and `v` to common size and add them.
"""
+ₚ(u::SVector{N1, T1}, v::SVector{N2, T2}) where {N1, T1, N2, T2} =
    vcat(u, @SVector zeros(T1, max(N1, N2) - N1)) +
    vcat(v, @SVector zeros(T2, max(N1, N2) - N2))

function +ₚ(u::AbstractVector, v::AbstractVector)
    n = length(u)
    m = length(v)
    vcat(u, zeros(T1, max(n, m) - n)) + vcat(v, zeros(T2, max(n, m) - m))
end

"""
    M *ₚ v

Promote matrix `M` and vector/matrix `v` to compatible sizes and multiply them.
TODO: better explanation.

# Examples

    julia> [0 0; 1 0; 0 1] *ₚ [1,1,1]
    4-element Array{Int64, 1}:
     0
     1
     1
     1
"""
function *ₚ(M::AbstractMatrix{T}, v::AbstractVector{U}) where {T, U}
    n = length(v)
    m1, m2 = size(M)
    if m2 ≥ n
        M * vcat(v, zeros(U, m2 - n))
    else
        d = n - m2
        [M zeros(T, m1, d);
         zeros(T, d, m2) I] * v
    end
end

function *ₚ(A::AbstractMatrix{T}, B::AbstractMatrix{U}) where {T, U}
    n1, n2 = size(A)
    m1, m2 = size(B)
    if m1 > n2
        d = m1 - n2
        [A zeros(T, n1, d);
         zeros(T, d, n2) I] * B
    elseif n2 > m1
        d = n2 - m1
        A * [B zeros(U, m1, d);
             zeros(U, d, m2) I]
    else
        A * B
    end
end
# = Abstract & helpers = #
#=
struct CurveProduct{D, T<:NTuple{D, Curve}} <: AbstractManifold{D}
    spaces ::T
end
Base.show(io::IO, cp::CurveProduct{D}) where D =
    print(io, join(string.(cp.spaces), " × "))

#CurveProduct(spaces::T) where {T<:NTuple{D, Curve}} where D =
#    CurveProduct{D, T}(spaces)

function offsetframe(ps::CurveProduct{D}, args::NTuple{D, T}) where {D, T}
    spaces = ps.spaces

    offset = zeros(T, D+3)
    frame  = Matrix{T}(I, (D+3, D+3))

    for i in 1:D
        space = spaces[i]
        scale = scaling(space)
        val, mat = offsetframe(space, args[i])
        offset += frame * vcat(val, zeros(T, D)) *
            (i > 1 ? scale(args[i-1]) : scale(zero(T)))
        frame *= idpad(mat, D+3)
    end
    offset, frame
end

(ps::CurveProduct{D})(args::Vararg{T, D}) where {D, T} = first(offsetframe(ps, args))

# All variants of constructing `CurveProduct`s.
LinearAlgebra.cross(c1::Curve, c2::Curve) =
    CurveProduct((c1, c2))
LinearAlgebra.cross(ps::CurveProduct{D}, c::Curve) where D =
    CurveProduct((ps.spaces..., c))
LinearAlgebra.cross(c::Curve, ps::CurveProduct{D}) where D =
    CurveProduct((c, ps.spaces))
LinearAlgebra.cross(p1::CurveProduct, p2::CurveProduct{D}) where D =
    CurveProduct((p1.spaces..., p2.spaces...))
=#

struct ProductSpace{D, C, T<:NTuple{2, AbstractManifold}} <: AbstractManifold{D, C}
    spaces::T
end
Base.show(io::IO, ps::ProductSpace) =
    print(io, "(", ps.spaces[1], " × ", ps.spaces[2], ")")

function LinearAlgebra.cross(m1::AbstractManifold{D1, C1},
                             m2::AbstractManifold{D2, C2}) where {D1,D2,C1,C2}
    D = D1 + D2
    C = max(D1 + C1, D1 + D2 + C2) - D
    spaces = (m1, m2)
    ProductSpace{D, C, typeof(spaces)}((m1, m2))
end

function offsetframe(ps::ProductSpace{D, C}, args::Vararg{T, D}) where {D, C, T}
    m1, m2 = ps.spaces

    val, mat = offsetframe(m1, args[1:dim(m1)]...)
    offset1 = vcat(val, zeros(T, D+C - length(val)))
    frame1  = idpad(mat, D+C)

    val, mat = offsetframe(m2, args[dim(m1).+(1:dim(m2))]...)
    offset2 = vcat(val, zeros(T, D+C - length(val)))
    frame2  = idpad(mat, D+C)

    offset1 + frame1 * offset2, frame1 * frame2
end

(ps::ProductSpace{D})(args::Vararg{T, D}) where {T, D} = offsetframe(ps, args...)[1]
