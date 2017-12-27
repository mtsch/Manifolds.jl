module Manifolds

using StaticArrays
using RecipesBase
using ForwardDiff: jacobian, Dual, partials, value

#include("transformations.jl")
#include("manifolds.jl")
include("frames.jl")
#include("plotting.jl")

#const Circle = NSphere{1}
#const Sphere = NSphere{2}

export
    # Transformations:
    basis, point, change_basis, translate, rotate_x, rotate_y, rotate_z,
    # Frame:
    tnbframe,
    # General:
    dim, codim, ambientdim,
    # Manifolds:
    Manifold,
    UnitSpace, Interval, Circle, Sphere, NSphere, Knot, ProductSpace, CartesianSpace
end
