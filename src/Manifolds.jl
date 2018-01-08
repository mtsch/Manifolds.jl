module Manifolds

using StaticArrays
using RecipesBase
using ForwardDiff: jacobian, Dual, partials, value

include("manifolds.jl")
include("frames.jl")
include("products.jl")
include("plotting.jl")

const Circle = NSphere{1}
const Sphere = NSphere{2}

export
    # Frame:
    tnbframe,
    # General:
    dim, codim, ambientdim,
    # Manifolds:
    Manifold,
    Interval, Circle, Sphere, NSphere, Knot,
    UnitSpace, ProductSpace, CartesianSpace
end
