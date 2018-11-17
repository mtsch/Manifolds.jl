using Manifolds
using Test

using LinearAlgebra
using Random
using StaticArrays
using Suppressor

# Randomized equality checking.
function ≃(m1::AbstractManifold{D}, m2::AbstractManifold{D}, n=1000) where D
    params = rand(D, n)
    vals1 = mapslices(x -> m1(x...), params, dims = 1)
    vals2 = mapslices(x -> m2(x...), params, dims = 1)
    if vals1 ≈ vals2
        true
    else
        bad = findall([vals1[:, i] ≉ vals2[:, i] for i in axes(vals1, 1)])
        println("Inconsistencies found in $(length(bad)) parameters. Examples:")
        foreach(Iterators.take(bad, 10)) do b
            println("  ", params[:, b], ", ", vals1[:, b], ", ", vals2[:, b])
        end
        false
    end
end

include("base.jl")
#include("modifiers.jl")
#include("products.jl")
