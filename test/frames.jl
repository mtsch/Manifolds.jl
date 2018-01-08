@testset "idpad" begin
    const idpad = Manifolds.idpad

    @test idpad([1, 2], 4) == [1, 2, 0, 0]
    @test idpad([1 2; 3 4], 3) == [1 2 0;
                                   3 4 0;
                                   0 0 1]
    @test idpad([1 2 3; 4 5 6], (3, 4)) == [1 2 3 0;
                                            4 5 6 0;
                                            0 0 0 1]
    @test idpad([1 2 3; 4 5 6], (4, 3)) == [1 2 3;
                                            4 5 6;
                                            0 0 0;
                                            0 0 0]

    sizes = 2.^(1:round(Int, sqrt(nruns)))
    for s1 in sizes, s2 in sizes
        v = rand(s1)
        m = rand(s1, s1)
        padv = idpad(v, s2)
        padm = idpad(m, s2)
        if s1 > s2
            @test padv == v[1:s2]
            @test padm == m[1:s2, 1:s2]
        elseif s1 < s2
            @test padv == [v; zeros(s2-s1)]
            @test padm == [m zeros(s1, s2-s1);
                           zeros(s2-s1, s1) I]
        else
            @test padv == v
            @test padm == m
        end

        @test size(padv) == (s2,)
        @test size(padm) == (s2, s2)

        @test IndexStyle(padv) == IndexStyle(typeof(padv)) == IndexLinear()
        @test IndexStyle(padm) == IndexStyle(typeof(padm)) == IndexCartesian()
    end
end

@testset "Frame" begin
    # dodaj test: tnbframe(Interval(), 0.) * [1,2,3,4,5]
    const Frame = Manifolds.Frame

    v3 = rand(3)
    v5 = rand(5)

    # Identities:
    @test Frame(eye(1), [0.0])      * v3 == v3
    @test Frame(eye(100), [0.0])    * v3 == v3
    @test Frame(eye(1), zeros(100)) * v3 == v3

    fr  = Frame(eye(3)[randcycle(3), :], zeros(3))
    id1 = Frame(eye(1), [0.0])
    id2 = Frame(eye(10), [0.0])
    id3 = Frame(eye(10), zeros(10))
    @test id1 * fr == fr * id1 == fr
    @test id2 * fr == fr * id2 == fr
    @test id3 * fr == fr * id3 == fr

    t = rand()
    tnbframe(Interval(), t) * v3 == [t; v3]
    tnbframe(Interval(), t) * v5 == [t; v5]
end
