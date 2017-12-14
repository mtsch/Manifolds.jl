module Manifolds

using StaticArrays
using RecipesBase
using HyperDualNumbers

include("transformations.jl")
include("curves.jl")
include("plotting.jl")

export
    # Transformations:
    change_basis, translate, scale, rotate_x, rotate_y, rotate_z,
    # Frame:
    frenetframe!, frenetframe,
    # Curves:
    circle, trefoil

end
