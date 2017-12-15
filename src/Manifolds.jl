module Manifolds

using StaticArrays
using RecipesBase
using HyperDualNumbers
using ForwardDiff: jacobian

include("transformations.jl")
include("manifolds.jl")
include("plotting.jl")

export
    # Transformations:
    basis, point, change_basis, translate, rotate_x, rotate_y, rotate_z,
    # Frame:
    frame,
    # General:
    dim, ambientdim,
    # Curves and surfaces:
    interval, circle, knot,
    torus, sphere, kleinbottle

end
