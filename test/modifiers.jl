@testset "modifiers.jl" begin
    @testset "reparametrized" begin
        interval = Cube(1)
        rp = reparametrized(interval, x->x^2)

        @test @capture_out(print(rp)) == @capture_err(print(stderr, rp))
        @test rp(0.0) == SVector(0.00)
        @test rp(0.5) == SVector(0.25)
        @test rp(1.0) == SVector(1.00)

        r = rand()
        @test Manifolds.offsetframe(rp, r) == Manifolds.offsetframe(interval, r^2)

        circle = Sphere(1)
    end

    @testset "scaled" begin
        interval = Cube(1)
        s = scaled(interval, 2)

        @test @capture_out(print(rm)) == @capture_err(print(stderr, rm))
        @test s(0.0) == SVector(0.0)
        @test s(0.5) == SVector(1.0)
        @test s(1.0) == SVector(2.0)

        @test 2interval == s

        r = rand()
        @test Manifolds.offsetframe(s, r) == Manifolds.offsetframe(interval, 2r)
    end

    @testset "translated" begin
    end

    @testset "combined" begin
    end
end


#= Reparametrized
show
call
construct

# Scaling
show
call
construct
*

# Translated
show
call
construct
+
=#
