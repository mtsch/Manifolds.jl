@testset "base.jl" begin
    @testset "Helpers and utilites" begin
        @testset "idpad" begin
            idpad = Manifolds.idpad

            a = rand(3, 2)
            @test idpad(a, 5) == [a zeros(3, 3); zeros(2, 2) I zeros(2)]
            a = rand(2, 3)
            @test idpad(a, 5) == [a zeros(2, 2); zeros(2, 3) I; zeros(1, 5)]
            @test idpad([1 0; 0 1], 10) == Matrix{Int}(I, 10, 10)

            @test_throws ArgumentError idpad(rand(10, 1), 5)
            @test_throws ArgumentError idpad(rand(1, 10), 5)
        end

        @testset "chopzeros" begin
            for i in 1:5
                a = [SVector{5, Float64}(vcat(rand(i), zeros(5-i))) for _ in 1:100]
                @test length(eltype(chopzeros(a))) == i
            end
        end
    end

    @testset "PointSpace" begin
        @test @capture_out(print(PointSpace())) == @capture_err(print(stderr, PointSpace()))

        @test iszero(PointSpace()())
        @test iszero(rand(PointSpace()))
        @test all(iszero, rand(PointSpace(), 100))

        @test Manifolds.hasoffsetframe(PointSpace())
        @test isempty(Manifolds.offsetframe(PointSpace())[1])
        @test isempty(Manifolds.offsetframe(PointSpace())[2])
    end

    @testset "ParametricCurve" begin
        circfun(x) = (cospi(2x), sinpi(2x))
        circ = ParametricCurve(circfun)

        @test @capture_out(print(circ)) == @capture_err(print(stderr, circ))

        @test norm(circ(rand())) ≈ 1
        @test norm(rand(circ)) ≈ 1
        @test all(norm.(rand(circ, 100)) .≈ 1)

        @test Manifolds.hasoffsetframe(circ)
        of = Manifolds.offsetframe(circ, 0.00)
        @test of[1] ≈ [1, 0, 0]
        @test of[2] ≈ [-1 0; 0 0; 0 1]
        of = Manifolds.offsetframe(circ, 0.25)
        @test of[1] ≈ [0, 1, 0]
        @test of[2] ≈ [0 0; -1 0; 0 1]
    end

    @testset "Sphere" begin
        for n in 1:5
            sphere = Sphere(n)
            @test @capture_out(print(sphere)) == @capture_err(print(stderr, sphere))

            @test norm(sphere(rand(n)...)) ≈ 1
            @test norm(rand(sphere)) ≈ 1
            @test all(norm.(rand(sphere, 100)) .≈ 1)

            if n == 1
                @test Manifolds.hasoffsetframe(sphere)
            else
                @test !Manifolds.hasoffsetframe(sphere)
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

            if n == 1
                @test Manifolds.hasoffsetframe(ball)
            else
                @test !Manifolds.hasoffsetframe(ball)
            end
        end
    end

    @testset "Cube" begin
        for n in 1:5
            cube = Cube(n)
            @test @capture_out(print(cube)) == @capture_err(print(stderr, cube))

            @test sum(cube(rand(n)...)) .≤ n
            @test sum(rand(cube)) ≤ n
            @test all(sum.(rand(cube, 100)) .≤ n)

            if n == 1
                @test Manifolds.hasoffsetframe(cube)
            else
                @test !Manifolds.hasoffsetframe(cube)
            end
        end
    end
end
