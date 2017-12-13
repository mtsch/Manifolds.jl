using StaticArrays

struct Transformation{D, T<:Real}
    basis::SMatrix{D, D, T}
    translate::SVector{D, T}
    scale::T
end

Transformation() = Transformation{1, Float64}(SMatrix{1,1,Float64}(1.0),
                                              SVector{1,Float64}(0.0),
                                              1.0)

function Transformation(m::AbstractMatrix, v::AbstractVector, s::S) where S
    # TODO: remove checks?
    size(m, 1) ≠ size(m, 2) &&
        throw(DimensionMismatch("`m` must be a square matrix!"))
    length(v) ≠ size(m, 1) &&
        throw(DimensionMismatch("`length(v)` must match the size of `m`!"))
    m'm ≠ eye(size(m, 1)) &&
        throw(ArgumentError("`m` must be orthogonal!"))


    T = promote_type(eltype(m), eltype(v), S)
    D = size(m, 1)
    Transformation{D, T}(SMatrix{D, D, T}(m),
                         SVector{D, T}(v),
                         convert(T, s))
end

function (t::Transformation{D})(arr::AbstractArray) where D
    size(arr, 2) ≠ 1 && size(arr, 1) ≠ size(arr, 2) &&
        throw(DimensionMismatch("`arr` must be a square matrix or vector!"))

    dim = max(D, size(arr, 1))

    m = idpad(t.basis, dim)
    v = idpad(t.translate, dim)
    a = idpad(arr, dim)

    (t.scale .* m) * a .+ v
end

function compose(t1::Transformation{D1, T1},
                 t2::Transformation{D2, T2}) where {D1, D2, T1, T2}

    T = promote_type(T1, T2)
    D = max(D1, D2)
    m1 = idpad(t1.basis, D)
    v1 = idpad(t1.translate, D)
    s1 = t1.scale
    m2 = idpad(t2.basis, D)
    v2 = idpad(t2.translate, D)
    s2 = t1.scale

    Transformation{D, T}(m1 * m2, t1(v2), s1 * s2)
end

Base.:∘(t1::Transformation, t2::Transformation) = compose(t1, t2)

function Base.inv(t::Transformation{D, T}) where {D, T}
    Transformation{D, T}(t.basis', -t.translate, 1/t.scale)
end

function Base.:(==)(t1::Transformation{D1},
                    t2::Transformation{D2}) where {D1, D2}
    D = max(D1, D2)
    m1 = idpad(t1.basis, D)
    v1 = idpad(t1.translate, D)
    s1 = t1.scale
    m2 = idpad(t2.basis, D)
    v2 = idpad(t2.translate, D)
    s2 = t1.scale

    m1 == m2 && v1 == v2 && s1 == s2
end

"""
    idpad(m, s::Int)

Pad a `SVector` to dimension `s` by extending them with zeros.
Pad a square `SMatrix` by expanding it with `I`.
"""
function idpad(m::SMatrix{D, D, T}, s::Int) where {D, T<:Real}
    #size(m, 1) ≠ size(m, 2) &&
    #    throw(DimensionMismatch("`m` must be a square matrix!"))

    Δ = s - D
    if Δ ≤ 0
        m
    else
        SMatrix{s, s}([m zeros(T, D, Δ);
                       zeros(T, Δ, D) I])
    end
end

function idpad(v::SVector{D, T}, s::Int) where {D, T<:Real}
    Δ = s - length(v)
    if Δ ≤ 0
        v
    else
        SVector{s}([v; zeros(T, Δ)])
    end
end

# ---------------------------------------------------------------------------- #
# Concrete transformations.
# ---------------------------------------------------------------------------- #
function change_basis(b::AbstractMatrix{T}) where T
    size(b, 1) ≠ size(b, 2) &&
        throw(DimensionMismatch("`basis` must be a square matrix!"))

    Transformation(b, zeros(SVector{size(b, 1), T}), one(T))
end

function translate(v::AbstractVector{T}) where T
    Transformation(eye(SMatrix{length(v), length(v), T}), v, one(T))
end

function scale(α::T) where T
    Transformation(eye(SMatrix{1,1,T}), zeros(SVector{1,T}), α)
end

rotate_x(θ::T) where T =
    change_basis(SMatrix{3,3,T}([1 0       0     ;
                                 0 cos(θ) -sin(θ);
                                 0 sin(θ)  cos(θ)]))
rotate_y(θ::T) where T =
    change_basis(SMatrix{3,3,T}([cos(θ) 0 sin(θ);
                                 0      1 0     ;
                                -sin(θ) 0 cos(θ)]))
rotate_z(θ::T) where T =
    change_basis(SMatrix{2,2,T}([cos(θ) -sin(θ);
                                 sin(θ)  cos(θ)]))

#TODO more rotations?
