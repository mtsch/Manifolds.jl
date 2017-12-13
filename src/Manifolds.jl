module Manifolds

using ForwardDiff: derivative
using StaticArrays
using RecipesBase

# package code goes here
include("transformations.jl")
include("curves.jl")


export
    # Frame:
    tangent, normal, binormal, frenetframe!, frenetframe,
    # Types:
    Circle, Trefoil

end # module
