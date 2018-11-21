abstract type ModifiedManifold{D, C} <: AbstractManifold{D, C} end

# reparametrize ========================================================================== #
struct ReparametrizedManifold{D, C, F, M<:AbstractManifold{D, C}} <: ModifiedManifold{D, C}
    base ::M
    par  ::F
end

Base.show(io::IO, rm::ReparametrizedManifold) =
    print(io, "reparametrized ($(rm.base))")

(sm::ReparametrizedManifold{D})(args::Vararg{T, D}) where {D, T} =
    sm.base(sm.par(args...)...)

reparametrized(man::AbstractManifold{D, C}, f) where {D, C} =
    ReparametrizedManifold{D, C, typeof(f), typeof(man)}(man, f)
reparametrized(man::AbstractManifold, t::Number) =
    reparametrized(man, x->x.*t)
reparametrized(::PointSpace, _) =
    PointSpace()

# scaling ================================================================================== #
# TODO? const scaling?
struct ScaledManifold{D, C, F, M<:AbstractManifold{D, C}} <: ModifiedManifold{D, C}
    base    ::M
    scaling ::F
end

Base.show(io::IO, sm::ScaledManifold) =
    print(io, "scaled ($(sm.base))")

scaling(::AbstractManifold) = one
scaling(mm::ModifiedManifold) = scaling(mm.base)
scaling(sm::ScaledManifold) = sm.scaling
scaling(sm::ScaledManifold{<:Any, <:Number}) = one

scaled(man::AbstractManifold{D, C}, f) where {D, C} =
    ScaledManifold{D, C, typeof(f), typeof(man)}(man, f)
scaled(::PointSpace, _) =
    PointSpace()
Base.:*(r::Number, man::AbstractManifold) =
    scaled(man, r)
Base.:*(man::AbstractManifold, r::Number) =
    scaled(man, r)

(sm::ScaledManifold{D, <:Any, <:Number})(args::Vararg{T, D}) where {T, D} =
    sm.scaling .* sm.base(args...)
(sm::ScaledManifold{D, <:Any})(args::Vararg{T, D}) where {T, D} =
    sm.base(args...)
Base.rand(rng::AbstractRNG, sm::ScaledManifold{D, <:Any, <:Number}) where D =
    sm.scaling .* rand(rng, sm.base)
Base.rand(rng::AbstractRNG, sm::ScaledManifold) =
    rand(rng, sm.base)

# translate ============================================================================== #
struct TranslatedManifold{D, C, T<:AbstractArray,
                          M<:AbstractManifold{D, C}} <: ModifiedManifold{D, C}
    base        ::M
    translation ::T
end

Base.show(io::IO, tm::TranslatedManifold) =
    print(io, "translated ($(tm.base))")

translation(::AbstractManifold) = false
translation(mm::ModifiedManifold) = translation(mm.base)
translation(tm::TranslatedManifold) = tm.translation

Base.:+(man::M, r::T) where {M<:AbstractManifold{D, C}, T<:AbstractArray} where {D, C} =
    TranslatedManifold{D, C, T, M}(man, r)
Base.:+(r, man::AbstractManifold{D}) where D = r + man

(tm::TranslatedManifold{D})(args::Vararg{T, D}) where {T, D} =
    promotedadd(tm.base(args...), tm.translation)
Base.rand(rng::AbstractRNG, tm::TranslatedManifold) =
    promotedadd(rand(rng, tm.base), tm.translation)
