using HyperDualNumbers
using StaticArrays
# General stuff
abstract type Curve{T<:Real} end

transformation(c::Curve, v) = c.transform * v .+ c.translate

function frenetframe!(res::AbstractArray, c::Curve, t)
    τ = hyper(t, 1.0, 1.0, 0.0)
    val = c(τ)
    res[:, 1] .= normalize!(eps1.(val))
    res[:, 2] .= eps1eps2.(val)
    res[:, 3] .= normalize!(res[:, 1] × res[:, 2])
    res[:, 2] .= -(res[:, 1] × res[:, 3])

    res
end

frenetframe(pc::Curve, t::T) where T =
    frenetframe!(Matrix{T}(3,3), pc, t)

Base.:*(m::SMatrix{3,3}, c::Curve) =
    T(m * c.transform, c.translate)
Base.:*(m::AbstractArray, c::Curve) =
    T(SMatrix{3, 3}(m) * c.transform, c.translate)
Base.:+(c::Curve, v::SVector) =
    T(c.transform, c.translate + v)
Base.:+(c::Curve, v::AbstractArray) =
    T(c.transform, c.translate + SVector{3}(v))

# Plotting
@recipe function f(c::Curve)
    pts = c.(linspace(0, 1, 1000))

    xs = getindex.(pts, 1)
    ys = getindex.(pts, 2)
    zs = getindex.(pts, 3)

    xs, ys, zs
end

struct CartesianCurves{N, T}
    curves::NTuple{N, Curve}
    transform::SMatrix{(3, 3), T}
    translate::SVector{3, T}
end

function Base.:×(args::Vararg{Curve{T}, N}) where {T, N}
    CartesianCurves{N, T}(convert(NTuple{N, Curve}, args),
                          eye(SMatrix{N+2, N+2}),
                          zeros(SVector{N+2}))
end

# ............................................................................ #
struct Circle{T} <: Curve{T}
    transform::SMatrix{(3, 3), T}
    translate::SVector{3, T}
end

# TODO boljši konstruktorji
Circle{T}() where T = Circle{T}(eye(SMatrix{3,3,T}), zeros(SVector{3,T}))
Circle(r::T) where T = Circle{T}(r*eye(SMatrix{3,3,T}), zeros(SVector{3,T}))
Circle() = Circle{Float64}()

(c::Circle)(t) = transformation(c, @SVector [cos(2π*t), sin(2π*t), zero(2π*t)])

struct Trefoil{T} <: Curve{T}
    transform::SMatrix{(3, 3), T}
    translate::SVector{3, T}
end

Trefoil{T}() where T = Trefoil(eye(SMatrx{3,3,T}), zeros(SVector{3,T}))
Trefoil() = Trefoil{Float64}()
(c::Trefoil)(t) = transformation(c, @SVector [sin(2π*t) + 2sin(4π*t),
                                              cos(2π*t) - 2cos(4π*t),
                                              -sin(6π*t)])

φ = π * rand()
rot = [cos(φ) 0 sin(φ);
       0      1 0
      -sin(φ) 0 cos(φ)]
trf = rot * Circle() + fill(2rand(), 3)
#trf = rot * Trefoil() + fill(2rand(), 3)
plt = plot(trf)
for t in linspace(0, 1, 5)
    point = trf[t]
    frame = frenetframe(trf, t)
    mframe = frame .+ point
    colors = [:red, :green, :blue]

    # Frenet frame.
    for i in 1:3
        plot!(plt,
              [point[1], mframe[1, i]],
              [point[2], mframe[2, i]],
              [point[3], mframe[3, i]],
              color=colors[i])
    end

    c = frame * Circle(0.5) + point
    plot!(c)
end
plt
