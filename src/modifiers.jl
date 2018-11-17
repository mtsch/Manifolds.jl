abstract type ModifiedManifold{D} <: AbstractManifold{D} end

# reparametrize ========================================================================== #
struct ReparametrizedManifold{D, F, M<:AbstractManifold} <: ModifiedManifold{D}
    base ::M
    par  ::F
end

Base.show(io::IO, rm::ReparametrizedManifold) =
    print(io, "reparametrized ($(rm.base))")

(sm::ReparametrizedManifold{D})(args::Vararg{T, D}) where {D, T} =
    sm.base(sm.par(args)...)

reparametrized(man::AbstractManifold{D}, f) where D =
    ReparametrizedManifold{D, typeof(f), typeof(man)}(man, f)
reparametrized(man::AbstractManifold, t::Number) =
    reparametrized(man, x->x.*t)
reparametrized(::PointSpace, _) =
    PointSpace()

# scaling ================================================================================== #
# TODO? const scaling?
struct ScaledManifold{D, F, M<:AbstractManifold} <: ModifiedManifold{D}
    base  ::M
    scaling ::F
end

Base.show(io::IO, sm::ScaledManifold) =
    print(io, "scaled ($(sm.base))")

scaling(::AbstractManifold) = one
scaling(mm::ModifiedManifold) = scaling(mm.base)
scaling(sm::ScaledManifold) = sm.scaling
scaling(sm::ScaledManifold{<:Any, <:Number}) = one

scaled(man::AbstractManifold{D}, f) where D =
    ScaledManifold{D, typeof(f), typeof(man)}(man, f)
scaled(::PointSpace, _) =
    PointSpace()
Base.:*(r::Number, man::AbstractManifold) =
    scaled(man, r)
Base.:*(man::AbstractManifold, r::Number) =
    scaled(man, r)

(sm::ScaledManifold{D, <:Number})(args::Vararg{T, D}) where {T, D} =
    sm.scaling .* sm.base(args...)
(sm::ScaledManifold{D})(args::Vararg{T, D}) where {T, D} =
    sm.base(args...)
Base.rand(rng::AbstractRNG, sm::ScaledManifold{D, <:Number}) where D =
    sm.scaling .* rand(rng, sm.base)
Base.rand(rng::AbstractRNG, sm::ScaledManifold) =
    rand(rng, sm.base)

# translate ============================================================================== #
struct TranslatedManifold{D, T<:AbstractArray, M<:AbstractManifold} <: ModifiedManifold{D}
    base        ::M
    translation ::T
end

Base.show(io::IO, tm::TranslatedManifold) =
    print(io, "translated ($(tm.base))")

translation(::AbstractManifold) = false
translation(mm::ModifiedManifold) = translation(mm.base)
translation(tm::TranslatedManifold) = tm.translation

Base.:+(man::M, r::T) where {M<:AbstractManifold{D}, T<:AbstractArray} where D =
    TranslatedManifold{D, T, M}(man, r)
Base.:+(r, man::AbstractManifold{D}) where D = r + man

function promotedadd(a::AbstractArray{T}, b::AbstractArray) where T
    da = length(a)
    db = length(b)
    if da > db
        vcat(SVector{db, T}(a[1:db] + b), SVector{da-db, T}(a[db+1:end]))
    else
        vcat(SVector{da, T}(a + b[1:da]), SVector{db-da, T}(b[da+1:end]))
    end
end

(tm::TranslatedManifold{D})(args::Vararg{T, D}) where {T, D} =
    promotedadd(tm.base(args...), tm.translation)
Base.rand(rng::AbstractRNG, tm::TranslatedManifold) =
    promotedadd(rand(rng, tm.base), tm.translation)
