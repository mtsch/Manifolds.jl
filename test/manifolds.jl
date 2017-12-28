# Number of points used in rand testing.
nrand = 100

@testset "Interval" begin
    @test dim(Interval(rand())) == 1
    @test codim(Interval(rand())) == 0
    @test ambientdim(Interval(rand())) == 1

    @testset "Basics" begin
        for _ in 1:nruns
            len = 5rand() + 1
            t   = rand()
            i1  = Interval(len)
            val = i1(t)
            @test length(val) == 1
            @test val[1] ≈ t * len

            i2  = Interval(x -> x^2)
            val = i2(t, scale = len)
            @test val[1] ≈ t * len^2
        end
    end

    @testset "rand" begin
        for _ in 1:nruns
            len = 5rand() + 1
            i1  = Interval(len)
            @test all(0 .< rand(i1, nrand) .< len)

            i2  = Interval(x -> x^2)
            @test all(0 .< rand(i1, nrand, scale = len) .< len^2)
        end
    end
end

@testset "NSphere" begin
    @testset "Basics" begin
        for d in 1:nruns
            rad = 5rand() + 1
            sp1 = NSphere{d}(rad)
            val = sp1(fill(rand(), d)...)
            @test length(val) == d + 1
            @test norm(val) ≈ rad

            sp2 = NSphere{d}(x -> x^2)
            val = sp2(fill(rand(), d)..., scale = rad)
            @test norm(val) ≈ rad^2

            @test dim(sp1) == d
            @test codim(sp1) == 1
            @test ambientdim(sp1) == d + 1
        end
    end

    @testset "rand" begin
        for d in 1:nruns
            rad = 5rand() + 1
            sp1 = NSphere{d}(rad)
            @test all(mapslices(norm, rand(sp1, nrand), 1) .≈ rad)

            sp2 = NSphere{d}(x -> x^2)
            @test all(mapslices(norm, rand(sp2, nrand, scale = rad), 1) .≈ rad^2)
        end
    end
end

@testset "Knot" begin
    # TODO: How to test this????
    @testset "Basics" begin
        for _ in 1:nruns
        end
    end
end

@testset "UnitSpace, copy" begin
    u = UnitSpace()
    for space in [Interval(),
                  Knot(),
                  NSphere{10}(),
                  Circle() * Circle(),
                  Circle() × Circle(),
                  Interval()^2 * (Knot() * NSphere{3}() × Interval() × Sphere()),
                  UnitSpace()]
        T = typeof(space)
        @test space == u * space == space * u == u * space * u
        @test space == u × space == space × u == u × space × u
        @test u == one(space) == one(T)

        @test copy(space) ≡ space
    end
end

@testset "ProductSpace" begin
    # Do Torus, Cylinder, Interval^n

    @testset "Torus" begin
        torus31 = Circle(3.0) * Circle()
        circ2   = Circle(2.0)
        circ3   = Circle(3.0)
        circ4   = Circle(4.0)

        for _ in 1:nruns
            t = rand()
            @test torus31(t, 0)    ≈ [circ2(t); 0] atol = 1e-16
            @test torus31(t, 0.5)  ≈ [circ4(t); 0] atol = 1e-16
            @test torus31(t, 0.25) ≈ [circ3(t); 1] atol = 1e-16
            @test torus31(t, 0.75) ≈ [circ3(t);-1] atol = 1e-16
        end
    end

    @testset "Intervals" begin
        for d in 1 + (1:nruns)
            lens = 5rand(d) + 1
            intd = prod(Interval(l) for l in lens)

            @test typeof(intd) == ProductSpace{d, 0, d}
            ts = rand(d)
            @test intd(ts...) ≈ ts .* lens
        end
    end

    @testset "Circles" begin
        for n in 1 + (1:nruns)
            mf = Circle()^n
            @test mf(fill(.25, n)...) ≈ vcat(0, fill(1, n))  atol = 1e-16
            @test mf(fill(.75, n)...) ≈ vcat(0, fill(-1, n)) atol = 1e-16
        end
    end

    @testset "tnbframe errors" begin
        # No error.
        circxglome = Circle() * NSphere{3}()
        @test circxglome            ≠ nothing
        @test Circle()^5            ≠ nothing
        @test Sphere()^5            ≠ nothing
        @test Interval()^5          ≠ nothing
        @test Circle() * circxglome ≠ nothing
        # Error / tnbframe not defined.
        @test_throws ErrorException NSphere{3}() * Circle()
        @test_throws ErrorException circxglome * Circle()
        @test_throws ErrorException circxglome * circxglome
    end
end

@testset "CartesianSpace" begin

    @testset "Clifford Torus" begin
        r1, r2 = 5rand(2) + 1
        clifft = Circle(r1) × Circle(r2)
        @test typeof(clifft) == CartesianSpace{2, 2, 2}

        for _ in 1:nruns
            t1, t2 = rand(2)
            @test clifft(t1, t2) ≈ [r1*cos(2π*t1), r1*sin(2π*t1),
                                    r2*cos(2π*t2), r2*sin(2π*t2)]
        end
    end

    @testset "Intervals" begin
        for d in 1 + (1:nruns)
            lens = 5rand(d) + 1
            intd = reduce(×, Interval(l) for l in lens)

            @test typeof(intd) == CartesianSpace{d, 0, d}
            ts = rand(d)
            @test intd(ts...) ≈ ts .* lens
        end
    end
end

@testset "when * == ×" begin
    for _ in 1:nruns
        t, u, v = rand(3)
        @test (Interval() * Interval())(t, u)  ≈ (Interval() × Interval())(t, u)
        @test (Interval() * Circle())(t, u)    ≈ (Interval() × Circle())(t, u)
        @test (Interval() * Sphere())(t, u, v) ≈ (Interval() × Sphere())(t, u, v)
    end
end
