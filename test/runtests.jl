using Manifolds
using Manifolds: +ₚ, *ₚ, offsetframe, dimincrease, dim, codim
using Test

using LinearAlgebra
using Random
using StaticArrays
using Suppressor

# Randomized equality checking. Also tests call and rand.
function ≅(m1::AbstractManifold{D}, m2::AbstractManifold{D}, n=1000) where D
    params = rand(D, n)
    vals1 = hcat(m1(zeros(D)...), mapslices(x -> m1(x...), params, dims = 1))
    vals2 = hcat(m2(zeros(D)...), mapslices(x -> m2(x...), params, dims = 1))
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
function ≅(m1::AbstractManifold, m2::AbstractManifold, n=1000)
    println("Mismatched dimensions.")
    false
end

⟂(v1::AbstractVector, v2::AbstractVector) = iszero(v1 ⋅ v2)

include("base.jl")
include("modifiers.jl")
#include("products.jl")
