@testset "Transformations" begin
    id = Manifolds.id(1)

    @testset "Rotations" begin
        for _ in 1:n_runs
            θ, φ = 2π * rand(2)

            @test rotate_x(θ) ∘ rotate_x(-θ) ≈ id
            @test rotate_y(θ) ∘ rotate_y(-θ) ≈ id
            @test rotate_z(θ) ∘ rotate_z(-θ) ≈ id

            @test rotate_x(θ) ∘ rotate_x(φ) ≈
                  rotate_x(φ) ∘ rotate_x(θ) ≈
                  rotate_x(θ + φ)
            @test rotate_y(θ) ∘ rotate_y(φ) ≈
                  rotate_y(φ) ∘ rotate_y(θ) ≈
                  rotate_y(θ + φ)
            @test rotate_z(θ) ∘ rotate_z(φ) ≈
                  rotate_z(φ) ∘ rotate_z(θ) ≈
                  rotate_z(θ + φ)

            @test rotate_x(θ) ∘ rotate_y(φ) ≠
                  rotate_y(φ) ∘ rotate_x(θ)
            @test rotate_x(θ) ∘ rotate_z(φ) ≠
                  rotate_z(φ) ∘ rotate_x(θ)
            @test rotate_y(θ) ∘ rotate_z(φ) ≠
                  rotate_z(φ) ∘ rotate_y(θ)

            @test rotate_x(θ)([1, 0, 0]) ≈ [1, 0, 0]
            @test rotate_y(θ)([0, 1, 0]) ≈ [0, 1, 0]
            @test rotate_z(θ)([0, 0, 1]) ≈ [0, 0, 1]
        end
        @test rotate_x(π/2)([0, 1 ,0]) ≈ [0, 0, 1]
        @test rotate_y(π/2)([0, 0, 1]) ≈ [1, 0, 0]
        @test rotate_z(π/2)([1, 0, 0]) ≈ [0, 1, 0]

        @test rotate_x(π/2)(eye(3)) ≈ [1  0  0;
                                       0  0 -1;
                                       0  1  0]
        @test rotate_y(π/2)(eye(3)) ≈ [0  0  1;
                                       0  1  0;
                                      -1  0  0]
        @test rotate_z(π/2)(eye(3)) ≈ [0 -1  0;
                                       1  0  0;
                                       0  0  1]

        @test rotate_x(2π) ≈ id
        @test rotate_y(2π) ≈ id
        @test rotate_z(2π) ≈ id
    end

    @testset "Translations" begin
        for dim in 1:n_runs
            a = 5rand(dim)
            b = 5rand(dim)

            @test translate(a) ∘ translate(-a) ≈
                  translate(-a) ∘ translate(a) ≈
                  id
            @test translate(a) ∘ translate(b) ≈
                  translate(b) ∘ translate(a) ≈
                  translate(a + b)

            @test translate(a)(b) ≈ a + b
            @test translate(b)(a) ≈ a + b
        end
    end

    @testset "Composite" begin
        tr = rotate_z(π/2) ∘  translate([0, 1, 0])
        @test tr([1, 0, 0]) ≈ [0, 2, 0]
        tr = translate([0,-1, 0]) ∘ tr
        @test tr([1, 0, 0]) ≈ [0, 1, 0]

        #TODO: more
    end
end
