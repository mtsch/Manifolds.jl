module Manifolds

using StaticArrays
using RecipesBase
using HyperDualNumbers

include("transformations.jl")
include("manifolds.jl")
include("plotting.jl")

export
    # Transformations:
    change_basis, translate, scale, rotate_x, rotate_y, rotate_z,
    # Frame:
    frame,
    # Curves:
    circle, trefoil

end
