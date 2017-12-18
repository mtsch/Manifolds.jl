idpad = Manifolds.idpad

@testset "ParametricCurves" begin
    @testset "sanity checks" begin
        @test dim(Interval()) == 1
        @test ambientdim(Interval()) == 1
        @test dim(Circle()) == 1
        @test ambientdim(Circle()) == 2
        @test dim(Knot()) == 1
        @test ambientdim(Knot()) == 3
    end

    curves = [interval(), circle(), knot()]

    @testset "frames" begin
        p, α = rand(2)
        f = frame(interval(α), p)
        @test basis(f) == eye(3)
        @test point(f) == [α*p, 0, 0]

        α = rand()
        f = frame(circle(α), 0.0)
        @test basis(f) ≈ [0 -1  0;
                          1  0  0;
                          0  0  1] atol = 1e-16
        @test point(f) ≈ [α, 0, 0] atol = 1e-16

        α = rand()
        f = frame(circle(α), 0.25)
        @test basis(f) ≈ [-1  0  0;
                           0 -1  0;
                           0  0  1] atol = 1e-16
        @test point(f) ≈ [ 0, α, 0] atol = 1e-16

        for t in rand(n_runs)
            f1 = frame(knot(), t)
            f2 = frame(knot(), t + 1)

            @test basis(f1) ≈ basis(f2) atol = 1e-16
            @test point(f1) ≈ point(f2) atol = 1e-16
        end

        for t in rand(n_runs)
            for c in curves
                @test idpad(c(t), 3) == idpad(point(frame(c, t)), 3)
            end
        end

    end

    @testset "eval" begin
        # Interval.
        for t in rand(n_runs)
            α = rand()
            @test interval(α)(t) == [α * t]
        end
        # Circle.
        for α in rand(n_runs)
            @test circle(α)(0)    ≈ [ α, 0]
            @test circle(α)(0.25) ≈ [ 0, α]
            @test circle(α)(0.5)  ≈ [-α, 0]
            @test circle(α)(0.75) ≈ [ 0,-α]
            @test circle(α)(1)    ≈ [ α, 0]
        end
        # Skipping knot.
    end
end
