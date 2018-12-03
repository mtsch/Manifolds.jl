# Frames ================================================================================= #
"""
    offsetframe(man::AbstractManifold{D, C}, args::NTuple{D, T})

Return a matrix that makes vectors perpendicular to `man(args...)` when multiplied with
`*ₚ` and the value of `man(args...)`.

Should return a `SVector{D+C}` and a matrix.
"""
offsetframe(man::AbstractManifold{D, C}, t::Vararg{T, D}) where {D, C, T} =
    man(t...), vcat(zeros(T, D+C, 1), one(T))
offsetframe(man::AbstractManifold{0}) =
    man(), zeros(0, 0)

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
            tangent × @SVector[0, 0, -1] : tangent × @SVector[0, -1, 0]
    end
    binormal = normalize(tangent × normal)
    normal   = normalize(binormal × tangent)

    offset[1:1+C], [normal binormal]
end

dimincrease(::AbstractManifold{0}) = 0
dimincrease(::AbstractManifold{D, C}) where {D, C} = D+C
dimincrease(::AbstractManifold{1, C}) where C = C ≤ 2 ? 1 : 1+C

"""
    u +ₚ v

Promote vectors `u` and `v` to common size and add them.
"""
+ₚ(u::SVector{N, T}, v::SVector{M, U}) where {N, T, M, U} =
    vcat(u, @SVector zeros(T, max(N, M) - N)) +
    vcat(v, @SVector zeros(U, max(N, M) - M))

function +ₚ(u::AbstractVector{T}, v::AbstractVector{U}) where {T, U}
    n = length(u)
    m = length(v)
    vcat(u, zeros(T, max(n, m) - n)) + vcat(v, zeros(U, max(n, m) - m))
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

# Product ================================================================================ #
struct ProductSpace{D, C, T<:NTuple{2, AbstractManifold}} <: AbstractManifold{D, C}
    spaces::T
end
Base.show(io::IO, ps::ProductSpace) =
    print(io, "(", ps.spaces[1], " × ", ps.spaces[2], ")")

# D, C calculation is wrong.
# TODO: add a trait dimensionincrease??
function LinearAlgebra.cross(m1::AbstractManifold{D1, C1},
                             m2::AbstractManifold{D2, C2}) where {D1,D2,C1,C2}
    D = D1 + D2
    #C = max(D1 + C1, D1 + D2 + C2) - D
    C = dimincrease(m1) - D1 + C2

    spaces = (m1, m2)
    ProductSpace{D, C, typeof(spaces)}((m1, m2))
end
LinearAlgebra.cross(::PointSpace, ::PointSpace) = PointSpace()
LinearAlgebra.cross(m::AbstractManifold, ::PointSpace) = m
LinearAlgebra.cross(::PointSpace, m::AbstractManifold) = m

function offsetframe(ps::ProductSpace{D, C}, args::Vararg{T, D}) where {D, C, T}
    m1, m2 = ps.spaces

    val1, mat1 = offsetframe(m1, args[1:dim(m1)]...)
    val2, mat2 = offsetframe(m2, args[dim(m1) .+ (1:dim(m2))]...)
    val2 *= scaling(m2)(args[dim(m1)]...)

    SVector{D+C, T}((val1 +ₚ mat1 *ₚ val2)[1:D+C]), mat1 *ₚ mat2
end

(ps::ProductSpace{D})(args::Vararg{T, D}) where {T, D} = offsetframe(ps, args...)[1]
