abstract type ModifiedManifold{D} <: AbstractManifold{D} end

struct ReparametrizedManifold{D, F, M<:AbstractManifold} <: ModifiedManifold{D}
    base ::M
    par  ::F
end

Base.show(io::IO, rm::ReparametrizedManifold) =
    print(io, "reparametrized $(rm.base)")

(sm::ReparametrizedManifold{D})(args::Vararg{T, D}) where {D, T} =
    sm.base(sm.par(args)...)

reparametrized(man::AbstractManifold{D}, f) where D =
    ReparametrizedManifold{D, typeof(f), typeof(man)}(man, f)
reparametrized(man::AbstractManifold, t::Real) =
    reparametrized(man, x->x.*t)

struct ScaledManifold{D, F, M<:AbstractManifold} <: ModifiedManifold{D}
    base  ::M
    scale ::F
end

Base.show(io::IO, sm::ScaledManifold) =
    print(io, "scaled $(sm.base)")

scale(::AbstractManifold) = _ -> 1
scale(mm::ModifiedManifold) = scale(mm.base)
scale(sm::ScaledManifold) = sm.scale
scale(sm::ScaledManifold{<:Any, <:Real}) = _ -> 1

scaled(man::AbstractManifold{D}, f) where D =
    ScaledManifold{D, typeof(f), typeof(man)}(man, f)
Base.:*(r::Real, man::AbstractManifold) =
    scaled(man, r)
Base.:*(man::AbstractManifold, r::Real) =
    scaled(man, r)

(sm::ScaledManifold{D, <:Real})(args::Vararg{T, D}) where {T, D} =
    sm.scale .* sm.base(args...)
(sm::ScaledManifold{D})(args::Vararg{T, D}) where {T, D} =
    sm.base(args...)
Base.rand(rng::AbstractRNG, sm::ScaledManifold{D, <:Real}, n::Integer) where D =
    sm.scale .* rand(rng, sm.base, n)
Base.rand(rng::AbstractRNG, sm::ScaledManifold, n::Integer) =
    rand(rng, sm.base, n)
