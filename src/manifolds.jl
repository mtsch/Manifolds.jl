# ============================================================================ #
abstract type Manifold{D, E} end

dim(::Manifold{D}) where D = D
ambientdim(::Manifold{D, E}) where {D, E} = E

# TODO: specialize for Manifold{1, 2}?
function frame(pc::Manifold{1}, t::T, tα=nothing) where T
    τ = hyper(t, one(T), one(T), zero(T))
    if tα == nothing
        pt = idpad(pc(τ), 3)
    else
        pt = idpad(pc(τ, tα), 3)
    end

    n = normalize(eps1.(pt))
    b = eps1eps2.(pt)
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), real.(pt))
end

function frame(pc::Manifold{1, 1}, t::T, tα=nothing) where T
    if tα == nothing
        RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t)[1],0,0]))
    else
        RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t, tα)[1],0,0]))
    end
end

function frame(ps::Manifold{2}, θ, φ, tα=nothing)
    if tα == nothing
        pt = idpad(ps(θ, φ), 3)
        J  = jacobian(v -> ps(v...), [θ, φ])
    else
        pt = idpad(ps(θ, φ, tα), 3)
        J  = jacobian(v -> ps(v...), [θ, φ, tα])
    end

    n = normalize(J[:, 1])
    b = J[:, 2]
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), pt)
end

# ============================================================================ #
struct ParametricCurve{E} <: Manifold{1, E}
    fun::Function
    scale::Function
end

# Constant scale
ParametricCurve(edim::Int, f::Function, scale = 1.0) =
    ParametricCurve{edim}(f, _ -> scale)

# Function scale
ParametricCurve(edim::Int, f::Function, scalefun::Function) =
    ParametricCurve{edim}(f, scalefun)

function (pc::ParametricCurve)(t, s = 0)
    t1, s1 = promote(t, s)
    pc.fun(t1) .* pc.scale(s1)
end

#=
# TODO: specialize for ParametricCurve{2}?
function frame(pc::ParametricCurve, t::T, ψ=0) where T
    τ = hyper(t, one(T), one(T), zero(T))
    #pt = pc.fun(τ) .* pc.scale
    pt  = pc(τ, ψ)
    pt3 = idpad(pt, 3)

    n = normalize(eps1.(pt3))
    b = eps1eps2.(pt3)
    t = normalize(n × b)
    b = -(n × t)

    RigidTransformation(SMatrix{3,3}([n b t]), real.(pt3))
end

function frame(pc::ParametricCurve{1}, t::T, α=0) where T
    RigidTransformation(eye(SMatrix{3,3,T}), @SVector([pc(t, α)[1],0,0]))
end
=#

# ============================================================================ #
struct ParametricSurface{E} <: Manifold{2, E}
    fun::Function
    scale::Function
end

# Constant scale
ParametricSurface(edim::Int, f::Function, scale = 1.0) =
    ParametricSurface{edim}(f, _ -> scale)

# Function scale
ParametricSurface(edim::Int, f::Function, scalefun::Function) =
    ParametricSurface{edim}(f, scalefun)

function (ps::ParametricSurface)(θ, φ, s = 0)
    t1, t2, s1 = promote(θ, φ, s)
    ps.fun(t1, t2) .* ps.scale(s1)
end

# ============================================================================ #
struct ProductSpace{D, E, N} <: Manifold{D, E}
    spaces::NTuple{N, Manifold}
end

#=
function gettrans(ps::ProductSpace{D, E, N},
                  args::Vararg{T, D}) where {T, D, E, N}
    spaces = ps.spaces

    argoffset = 0
    trans = id(1, T)
    α = zero(T)
    for i in 1:N-1
        d = dim(spaces[i])
        trans = lift(trans(frame(spaces[i], args[argoffset + (1:d)]..., α)), d)

        α = args[argoffset + 1]
        argoffset += d
    end
    trans, argoffset, α
end
=#

(ps::ProductSpace{D})(args::Vararg{Any, D}) where D =
    ps(promote(args...)...)
function (ps::ProductSpace{D, E, N})(args::Vararg{T, D}) where {T, D, E, N}
    #trans, argoffset, α = gettrans(ps, args...)
    spaces = ps.spaces

    argoffset = 0
    trans = id(1, T)
    t = zero(T)
    for i in 1:N-1
        d = dim(spaces[i])
        trans = lift(trans(frame(spaces[i], args[argoffset + (1:d)]..., t)), d)

        t = args[argoffset + 1]
        argoffset += d
    end

    trans(ps.spaces[end](args[argoffset+1:end]..., t))[1:E]
end

#=
frame(ps::ProductSpace{D}, args::Vararg{Any, D}) where D =
    frame(ps, promote(args...)...)
function frame(ps::ProductSpace{D, E, N}, args::Vararg{T, D}) where {T, D, E, N}
    trans, argoffset, α = gettrans(ps, args...)
    trans(frame(ps.spaces[end], args[argoffset+1:end]..., α))
end
=#

function Base.:*(m1::Manifold{D1, E1},
                 m2::Manifold{D2, E2}) where {D1, D2, E1, E2}
    ProductSpace{D1+D2, max(E1, D1+E2), 2}((m1, m2))
end

function Base.:*(ps::ProductSpace{D1, E1, N},
                 m::Manifold{D2, E2}) where {D1, D2, E1, E2, N}
    ProductSpace{D1+D2, max(E1, D1+E2), N+1}((ps.spaces..., m))
end

function Base.:*(m::Manifold{D1, E1},
                 ps::ProductSpace{D2, E2, N}) where {D1, D2, E1, E2, N}
    ProductSpace{D1+D2, max(E1, D1+E2), N+1}((m, ps.spaces...))
end

function Base.:*(ps1::ProductSpace{D1, E1, N1},
                 ps2::ProductSpace{D2, E2, N2}) where {D1, D2, E1, E2, N1, N2}
    ProductSpace{D1+D2, max(E1, D1+E2), N1+N2}((ps1.spaces..., ps1.spaces...))
end


# ============================================================================ #
# == EXAMPLES == #
interval(α=1.0) =
    ParametricCurve(1, t -> SVector(t), α)

circle(r=1.0) =
    ParametricCurve(2, t -> @SVector([cos(2π*t), sin(2π*t)]), r)

knot(p=2, q=3, n=1, m=1.5, h=1, α=1) =
    ParametricCurve(3, t -> @SVector([m*cos(2π*p*t) + n*cos((2-q)*2π*t),
                                      m*sin(2π*p*t) + n*sin((2-q)*2π*t),
                                      h * sin(q*2π*t)]), α)

torus2(R, r, α=1.0) =
    ParametricSurface(3, (θ, φ) -> @SVector([(R + r*cos(2π*θ))*cos(2π*φ),
                                             (R + r*cos(2π*θ))*sin(2π*φ),
                                             r*sin(2π*θ)]), α)

sphere(α=1.0) =
    ParametricSurface(3, (θ, φ) -> @SVector([cos(2π*φ)*sin(π*θ),
                                             sin(2π*φ)*sin(π*θ),
                                             cos(π*θ)]), α)

kleinbottle(R, P, ϵ=0, α=1.0) =
    ParametricSurface(4, (θ, φ) ->
                      @SVector([R*(cos(π*θ)*cos(2π*φ) - sin(π*θ)*sin(4π*φ)),
                                R*(sin(π*θ)*cos(2π*φ) - cos(π*θ)*sin(4π*φ)),
                                P*cos(2π*θ) * (1 + ϵ*sin(2π*φ)),
                                P*sin(2π*θ) * (1 + ϵ*sin(2π*φ))]), α)


# == GARBAGE == #
if false

    #torus = trefoil() * circle(.2)
    surf = torus2(3.,1.)#circle(3.) * circle()#torus2(1.,.1)#torus2(1., .1)
    plt = plot(surf, alpha=0.5)
    for i in 0:4
        #θ, φ = rand(2)#0., 1.
        θ = 0
        φ = i / 4
        p = surf(θ, φ)
        plot!(plt, frame(surf, θ, φ))
        @show basis(frame(surf, θ, φ))
    end
    plt

end
