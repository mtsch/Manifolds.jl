# TODO: * generalize products:
#           * everything uses offsetframes
#           * add flat product -- ModifiedManifold that prevents offsetframes from forming
#           * add default flat offset frame
#       * AbstractManifold{D, C} -- codimension
#       * ParametricCurve --> ParametricManifold{D, C}, curves are special cases
module Manifolds

using LinearAlgebra
using Random
using ForwardDiff: Dual, partials, value
using RecipesBase
using StaticArrays

include("base.jl")
include("modifiers.jl")
include("products.jl")
include("plotting.jl")

export
    # Types
    AbstractManifold,
    # Primitives
    PointSpace, ParametricCurve, Sphere, Ball,
    # Modifiers
    ModifiedManifold, ReparametrizedManifold, ScaledManifold, TranslatedManifold,
    reparametrized, scaled, translated,
    # Helpers
    chopzeros
end
