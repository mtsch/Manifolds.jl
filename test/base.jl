@testset "base.jl" begin
    @testset "PointSpace" begin
        @test @capture_out(print(PointSpace())) == @capture_err(print(stderr, PointSpace()))
        @test iszero(PointSpace()())
        @test iszero(rand(PointSpace()))
        @test all(iszero, rand(PointSpace(), 100))
        @test dim(PointSpace()) == 0
        @test codim(PointSpace()) == 1

        @test offsetframe(PointSpace())[1] == PointSpace()()
        @test isempty(offsetframe(PointSpace())[2])
        r = rand(5)
        @test offsetframe(PointSpace())[2] *ₚ r == r
        @test dimincrease(PointSpace()) == 0
    end

    @testset "Sphere" begin
        for n in 1:5
            sphere = Sphere(n)
            @test @capture_out(print(sphere)) == @capture_err(print(stderr, sphere))
            @test norm(sphere(rand(n)...)) ≈ 1
            @test norm(rand(sphere)) ≈ 1
            @test all(norm.(rand(sphere, 100)) .≈ 1)
            @test dim(sphere) == n
            @test codim(sphere) == 1
            if n == 1
                @test dimincrease(sphere) == 1
                of = offsetframe(sphere, 0.00)
                @test of[1] ≈ [1, 0]
                @test of[2] ≈ [-1 0; 0 0; 0 1]
                of = offsetframe(sphere, 0.25)
                @test of[1] ≈ [0, 1]
                @test of[2] ≈ [0 0; -1 0; 0 1]
            else
                @test dimincrease(sphere) == n+1
                _, f = offsetframe(sphere, rand(n)...)
                r = rand(5)
                @test f *ₚ r == vcat(zeros(n+1), r)
            end
        end
    end

    @testset "Ball" begin
        for n in 1:5
            ball = Ball(n)
            @test @capture_out(print(ball)) == @capture_err(print(stderr, ball))
            @test norm(ball(rand(n)...)) ≤ 1
            @test norm(rand(ball)) ≤ 1
            @test all(norm.(rand(ball, 100)) .≤ 1)
            @test dim(ball) == n
            @test codim(ball) == 0
            @test dimincrease(ball) == n
            r = rand(5)
            @test offsetframe(ball, rand(n)...)[2] *ₚ r == vcat(zeros(n), r)
        end
    end

    @testset "Cube" begin
        for n in 1:5
            cube = Cube(n)
            @test @capture_out(print(cube)) == @capture_err(print(stderr, cube))
            @test sum(cube(rand(n)...)) .≤ n
            @test sum(rand(cube)) ≤ n
            @test all(sum.(rand(cube, 100)) .≤ n)
            @test dim(cube) == n
            @test codim(cube) == 0
            @test dimincrease(cube) == n
            r = rand(5)
            @test offsetframe(cube, rand(n)...)[2] *ₚ r == vcat(zeros(n), r)
        end
    end

    @testset "ParametricManifold" begin
        circfun(x) = (cospi(2x), sinpi(2x))
        circ = ParametricManifold{1, 1}(circfun)
        @test @capture_out(print(circ)) == @capture_err(print(stderr, circ))
        @test norm(circ(rand())) ≈ 1
        @test norm(rand(circ)) ≈ 1
        @test all(norm.(rand(circ, 100)) .≈ 1)
        @test dim(circ) == 1
        @test codim(circ) == 1
        of = offsetframe(circ, 0.00)
        @test of[1] ≈ [1, 0]
        @test of[2] ≈ [-1 0; 0 0; 0 1]
        of = offsetframe(circ, 0.25)
        @test of[1] ≈ [0, 1]
        @test of[2] ≈ [0 0; -1 0; 0 1]
        @test dimincrease(circ) == 1
        r = rand(5)
        @test length(offsetframe(circ, rand())[2] *ₚ r) == 6
        @test circ ≅ Sphere(1)

        # TODO: others
    end
end
