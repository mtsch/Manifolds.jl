using Manifolds
using Base.Test

curves = [Circle(), Trefoil()]

@testset "Frenet Frame" begin
    for c in curves, _ in 1:10
        t = rand()

        frame = frenetframe(c, t)
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
    for φ in [0; rand(9)]
        rot = [cospi(φ) 0 sinpi(φ);
               0        1 0
              -sinpi(φ) 0 cospi(φ)]

        c = rot * Circle() + fill(rand(), 3)

        frame = frenetframe(c, 0.0)
        @test frame ≈ rot * [0 -1  0;
                             1  0  0;
                             0  0  1] atol = 1e-16

        frenetframe!(frame, c, 0.25)
        @test frame ≈ rot * [-1  0  0;
                              0 -1  0;
                              0  0  1] atol = 1e-16

        frenetframe!(frame, c, 0.5)
        @test frame ≈ rot * [0  1  0;
                            -1  0  0;
                             0  0  1] atol = 1e-16

        frenetframe!(frame, c, 0.75)
        @test frame ≈ rot * [1  0  0;
                             0  1  0;
                             0  0  1] atol = 1e-16

        frenetframe!(frame, c, 1.0)
        @test frame ≈ rot * [0 -1  0;
                             1  0  0;
                             0  0  1] atol = 1e-16
    end

end
