@testset "modifiers.jl" begin
    @testset "reparametrized" begin
        # Interval by x->x²
        interval = Cube(1)
        ri = reparametrized(interval, x->x^2)
        @test @capture_out(print(ri)) == @capture_err(print(stderr, ri))
        @test ri(0.0) == SVector(0.00)
        @test ri(0.5) == SVector(0.25)
        @test ri(1.0) == SVector(1.00)
        r = rand()
        @test offsetframe(ri, r)[1] ≈ offsetframe(interval, r^2)[1]
        @test offsetframe(ri, r)[2] ≈ offsetframe(interval, r^2)[2]

        # Sphere by 2
        circle = Sphere(1)
        rc = reparametrized(circle, 2)
        @test @capture_out(print(rc)) == @capture_err(print(stderr, rc))
        @test rc(0.0) == rc(0.5) == rc(1.0)
        r = rand()
        @test offsetframe(rc, r)[1] ≈ offsetframe(circle, 2r)[1]
        @test offsetframe(rc, r)[2] ≈ offsetframe(circle, 2r)[2]
        @test norm(rand(rc)) ≈ 1

        # Identity reprametrization does not change the space.
        for i in 1:5
            @test reparametrized(Sphere(i), (args...) -> args) ≅ Sphere(i)
        end
    end

    @testset "scaled" begin
        #
        interval = Cube(1)
        si = scaled(interval, 2)
        @test @capture_out(print(si)) == @capture_err(print(stderr, si))
        @test si(0.0) == SVector(0.0)
        @test si(0.5) == SVector(1.0)
        @test si(1.0) == SVector(2.0)
        @test 2interval == si
        r = rand()
        @test offsetframe(si, r) == offsetframe(interval, 2r)
        @test scaled(interval, 2) == 2interval == interval * 2

        # Radius of sphere should match scaling.
        for i in 1:5
            r = rand()
            sc = r * Sphere(i)
            @test norm(rand(sc)) ≈ r
        end

        # Mulitplying by 1 doesn't change the space.
        for i in 1:5
            @test 1 * Sphere(i) ≅ Sphere(i)
        end
    end

    @testset "translated" begin
    end

    @testset "flattened" begin
        fs = flattened(Sphere(1))
        @test @capture_out(print(fs)) == @capture_err(print(stderr, fs))
        @test offsetframe(fs, 0.0)[1] ≈ SVector(1.0, 0.0)
        @test offsetframe(fs, 0.0)[2] ≈ reshape([0.0, 0.0, 1.0], (3, 1))
        @test fs ≅ Sphere(1)

        # Noop for spaces higher dimensions.
        for i in 2:5
            fi = flattened(Sphere(i))
            r = rand(i)
            @test offsetframe(fi, r...)[1] ≈ offsetframe(Sphere(i), r...)[1]
            @test offsetframe(fi, r...)[2] ≈ offsetframe(Sphere(i), r...)[2]
            @test fi ≅ Sphere(i)
        end
    end
end
