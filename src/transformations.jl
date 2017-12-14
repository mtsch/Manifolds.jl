# TODO: Remove scale, rename to RigidRigidTransform

struct RigidTransform{D, T<:Real}
    basis::SMatrix{D, D, T}
    translate::SVector{D, T}
end

RigidTransform() = RigidTransform{1, Float64}(SMatrix{1,1,Float64}(1.0),
                                              SVector{1,Float64}(0.0))

function RigidTransform(m::AbstractMatrix, v::AbstractVector)
    # TODO: remove checks?
    size(m, 1) ≠ size(m, 2) &&
        throw(DimensionMismatch("`m` must be a square matrix!"))
    length(v) ≠ size(m, 1) &&
        throw(DimensionMismatch("`length(v)` must match the size of `m`!"))
    m'm ≠ eye(size(m, 1)) &&
        throw(ArgumentError("`m` must be orthogonal!"))


    T = promote_type(eltype(m), eltype(v))
    D = size(m, 1)
    RigidTransform{D, T}(SMatrix{D, D, T}(m),
                         SVector{D, T}(v))
end

function (t::RigidTransform{D})(arr::AbstractArray, translate=true) where D
    # TODO: remove checks?
    size(arr, 2) ≠ 1 && size(arr, 1) ≠ size(arr, 2) &&
        throw(DimensionMismatch("`arr` must be a square matrix or vector!"))

    dim = max(D, size(arr, 1))

    m = idpad(t.basis, dim)
    a = idpad(arr, dim)

    if translate
        v = idpad(t.translate, dim)
        m * a .+ v
    else
        m * a
    end
end

function compose(t1::RigidTransform{D1, T1},
                 t2::RigidTransform{D2, T2}) where {D1, D2, T1, T2}

    T = promote_type(T1, T2)
    D = max(D1, D2)
    m1 = idpad(t1.basis, D)
    v1 = idpad(t1.translate, D)
    m2 = idpad(t2.basis, D)
    v2 = idpad(t2.translate, D)

    RigidTransform{D, T}(m1 * m2, v1 + v2)
end

Base.:∘(t1::RigidTransform, t2::RigidTransform) = compose(t1, t2)

function Base.inv(t::RigidTransform{D, T}) where {D, T}
    RigidTransform{D, T}(t.basis', -t.translate)
end

function Base.:(==)(t1::RigidTransform{D1},
                    t2::RigidTransform{D2}) where {D1, D2}
    D = max(D1, D2)
    m1 = idpad(t1.basis, D)
    v1 = idpad(t1.translate, D)
    m2 = idpad(t2.basis, D)
    v2 = idpad(t2.translate, D)

    m1 == m2 && v1 == v2
end

function Base.isapprox(t1::RigidTransform{D1},
                       t2::RigidTransform{D2};
                       kwargs...) where {D1, D2}
    D = max(D1, D2)
    m1 = idpad(t1.basis, D)
    v1 = idpad(t1.translate, D)
    m2 = idpad(t2.basis, D)
    v2 = idpad(t2.translate, D)

    isapprox(m1, m2; kwargs...) && isapprox(v1, v2; kwargs...)
end

"""
    idpad(m, s::Int)

Pad a `SVector` to dimension `s` by extending them with zeros.
Pad a square `SMatrix` by expanding it with `I`.
"""
function idpad(m::SMatrix{D, D, T}, s::Int) where {D, T}

    Δ = s - D
    if Δ ≤ 0
        m
    else
        SMatrix{s, s}([m zeros(T, D, Δ);
                       zeros(T, Δ, D) I])
    end
end

function idpad(v::SVector{D, T}, s::Int) where {D, T}
    Δ = s - length(v)
    if Δ ≤ 0
        v
    else
        SVector{s}([v; zeros(T, Δ)])
    end
end

function idpad(m::AbstractMatrix, s::Int)
    size(m, 1) ≠ size(m, 2) &&
        throw(DimensionMismatch("`m` must be a square matrix!"))
    idpad(SMatrix{size(m, 1), size(m, 1)}(m), s)
end

function idpad(v::AbstractVector, s::Int)
    idpad(SVector{length(v)}(v), s)
end


# ---------------------------------------------------------------------------- #
# Concrete transformations.
# ---------------------------------------------------------------------------- #
function change_basis(b::AbstractMatrix{T}) where T
    size(b, 1) ≠ size(b, 2) &&
        throw(DimensionMismatch("`basis` must be a square matrix!"))

    RigidTransform(b, zeros(SVector{size(b, 1), T}))
end

function translate(v::AbstractVector{T}) where T
    RigidTransform(eye(SMatrix{length(v), length(v), T}), v)
end

id(d::Int, T::Type=Int) = RigidTransform{d, T}(eye(SMatrix{d,d,T}),
                                               zeros(SVector{d,T}))

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

#TODO more rotations/general rotation?
