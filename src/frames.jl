struct IDPadding{T, N, A<:AbstractArray{T, N}} <: AbstractArray{T, N}
    data::A
    size::NTuple{N, Int}
end

Base.size(idp::IDPadding) = idp.size

Base.IndexStyle(::IDPadding{T, 1}) where T = IndexLinear()
Base.IndexStyle(::IDPadding{T, 2}) where T = IndexCartesian()

function Base.getindex(idp::IDPadding{T, 1}, i) where T
    if i > length(idp.data)
        zero(T)
    else
        idp.data[i]
    end
end

function Base.getindex(idp::IDPadding{T, 2}, i, j) where T
    s1, s2 = size(idp.data)

    if i ≤ s1 && j ≤ s2
        idp.data[i, j]
    elseif i - s1 == j - s2
        one(T)
    else
        zero(T)
    end
end

"""
    idpad(arr::AbstractArray{T, N}, size::NTuple{N, Int})

Pad array `arr` with idenity matrix of the same dimension to `size`.
Pad vectors with zeros.
"""
idpad(arr::AbstractMatrix{T}, size::Tuple{Int, Int}) where T =
    IDPadding{T, 2, typeof(arr)}(arr, size)
idpad(arr::AbstractMatrix{T}, size::Int) where T =
    IDPadding{T, 2, typeof(arr)}(arr, (size, size))
idpad(arr::AbstractVector{T}, size::Int) where T =
    IDPadding{T, 1, typeof(arr)}(arr, (size,))

# ============================================================================ #
# TODO: clean up the chompzeros mess and add inplace Frame composition.
"""
    Frame(basis_change, translation)

The frame object represents a change of basis and a translation. It can be
multiplied with a vector much like a matrix to transform the vector. The frame
can be applied to a vector of any length.
"""
struct Frame{T, M<:AbstractMatrix{T}, V<:AbstractVector{T}}
    basis_change::M
    translation::V
end

function chompzeros!(vec::AbstractVector)
    newlength = length(vec)
    for i in length(vec):-1:2
        vec[i] ≠ 0 && break
        newlength -= 1
    end
    resize!(vec, newlength)
end

function chompzeros!(vec::SVector)
    newlength = length(vec)
    for i in length(vec):-1:2
        vec[i] ≠ 0 && break
        newlength -= 1
    end
    # TODO always force SVector?
    vec[1:newlength]
end

function chompzeros(mat::AbstractMatrix)
    new_h = size(mat, 1)
    new_w = size(mat, 2)
    while new_h > 1 && all(iszero.(mat[new_h, :]))
        new_h -= 1
    end
    while new_w > 1 && all(iszero.(mat[:, new_w]))
        new_w -= 1
    end
    mat[1:new_h, 1:new_w]
end

function Frame(basis::AbstractMatrix{T},
               translation::AbstractVector{T}) where T
    chomped = chompzeros!(translation)
    M = typeof(basis)
    V = typeof(chomped)
    Frame{T, M, V}(basis, chomped)
end

function Base.:(==)(fr1::Frame, fr2::Frame)
    size1 = size(fr1.basis_change)
    size2 = size(fr2.basis_change)
    common_size = (max(size1[1], size2[1]), max(size1[2], size2[2]))

    idpad(fr1.basis_change, common_size) == idpad(fr2.basis_change, common_size) &&
    fr1.translation == fr2.translation
end

function Base.:*(fr::Frame{T1}, vec::AbstractVector{T2}) where {T1, T2}
    dim_incr  = size(fr.basis_change, 1) - size(fr.basis_change, 2)
    in_dim    = length(vec)
    out_dim   = in_dim + dim_incr
    trans_dim = length(fr.translation)

    T = promote_type(T1, T2)
    res = zeros(T, max(out_dim, trans_dim))
    res[1:out_dim] .= idpad(fr.basis_change, (out_dim, in_dim)) * vec
    res[1:trans_dim] .+= fr.translation
    res
end

function Base.:*(fr1::Frame{T1}, fr2::Frame{T2}) where {T1, T2}
    T = promote_type(T1, T2)

    out1, in1 = size(fr1.basis_change)
    out2, in2 = size(fr2.basis_change)
    common_dim = max(in1, in2, out1, out2)
    padding = (common_dim + 1, common_dim + 1)

    new_basis = idpad(fr1.basis_change, padding) *
        idpad(fr2.basis_change, padding)

    Frame(chompzeros(new_basis), fr1 * fr2.translation)
end

# ============================================================================ #
"""
    tnbframe(man::Manifold, ts...; scale=nothing)

Return a coordinate system (wrapped in a `Frame` object) that is perpendicular
to the manifold `man` evaluated at `(ts, scale)`.

For example for a curve, the `Frame` is a mapping:

* x ↦ N
* y ↦ B

where N and B are the normal and binormal vectors.
z and higher coordinates are mapped to unit vectors in higher dimensions, in the
curve case, z ↦ [0,0,0,1].

The frame also includes a translation that is equal to
`man(ts..., scale = scale)`.
"""
function tnbframe(man::Manifold{1, C}, t::T, scale=nothing) where {C, T}
    C > 2 && error("Cannot calculate the tnbframe!") # TODO: figure this out - it might work.

    t_dual = Dual{:d1}(Dual{:d2}(t, one(T)), one(T))
    if scale == nothing
        val = collect(idpad(man(t_dual), 3))
    else
        val = collect(idpad(man(t_dual, scale=scale), 3))
    end

    tangent  = normalize(value.(partials.(val, 1)))
    normal   = partials.(partials.(val, 1), 1)
    binormal = normalize(tangent × normal)
    normal   = binormal × tangent

    Frame([normal binormal], value.(value.(val)))
end

function tnbframe(man::Manifold{1, 0}, t, scale=nothing)
    if scale == nothing
        val = man(t)
    else
        val = man(t, scale = scale)
    end
    T = eltype(val)
    Frame(reshape([zero(T), one(T)], 2, 1), val)
end

function tnbframe(man::Manifold{2}, t1, t2, scale=nothing)
    if scale == nothing
        val = idpad(man(t1, t2), 3)
        J   = jacobian(v -> man(v...), [t1, t2])
    else
        val = idpad(man(t1, t2, scale = scale), 3)
        J   = jacobian(v -> man(v[1:end-1]..., scale=v[end]), [t1, t2, scale])
    end

    normal = normalize(J[:, 2] × J[:, 1])

    Frame(reshape(normal, length(binormal), 1), val)
end
