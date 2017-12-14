@testset "Curves" begin
    curves = [circle(1.0), trefoil(1.0)]

    @testset "Frenet Frame" begin
        for c in curves, θ in [0; rand(n_runs-1)]
            t = rand()
            rot = rotate_y(2π * θ)

            frame, _ = frenetframe(c, t)
            N = frame[:, 1]
            T = frame[:, 2]
            B = frame[:, 3]

            # Frame vectors are normal.
            @test norm(N) ≈ 1
            @test norm(T) ≈ 1
            @test norm(B) ≈ 1

            # Frame vectors are orthogonal.
            # Testing for a ≈ 0 does not work.
            @test N ⋅ T + 1 ≈ 1
            @test N ⋅ B + 1 ≈ 1
            @test B ⋅ T + 1 ≈ 1
        end

        # Check directions on a rotated circle.
        for θ in [0; rand(n_runs-1)]
            rot = rotate_y(2π * θ)

            c = rot(circle(1.0))

            frame, _ = frenetframe(c, 0.0)
            @test frame ≈ rot([0 -1  0;
                               1  0  0;
                               0  0  1]) atol = 1e-16

            frame, _ = frenetframe(c, 0.25)
            @test frame ≈ rot([-1  0  0;
                                0 -1  0;
                                0  0  1]) atol = 1e-16

            frame, _ = frenetframe(c, 0.5)
            @test frame ≈ rot([0  1  0;
                              -1  0  0;
                               0  0  1]) atol = 1e-16

            frame, _ = frenetframe(c, 0.75)
            @test frame ≈ rot([1  0  0;
                               0  1  0;
                               0  0  1]) atol = 1e-16

            frame, _ = frenetframe(c, 1.0)
            @test frame ≈ rot([0 -1  0;
                               1  0  0;
                               0  0  1]) atol = 1e-16
        end
    end
end
