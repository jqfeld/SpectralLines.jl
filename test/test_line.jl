using Test, SpectralLines

using QuadGK

const test_points = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 100., 1000., 1e5]
const test_shapes = [Gaussian(1.0), Lorentzian(1.0), Voigt(1.0, 0.5), PseudoVoigt(1.0, 0.5)]

@testset "Gaussian Line Shape" begin
    @testset "Amplitude factor" begin
        test_amplitudes = [1e-6, 1.0, 1e6, 1e10]
        for shape in test_shapes
            for ampl in test_amplitudes
                l = Line(0., ampl, shape)
                for x in test_points
                    @test l(x) ≈ ampl * shape(x) rtol = 1e-14
                end
            end
        end
    end

    @testset "Position" begin
        test_positions = [0.1, 1., 10, 100, 500., 1e4]
        for shape in test_shapes
            for pos in test_positions
                lp = Line(pos, 1., shape)
                ln = Line(-pos, 1., shape)
                @test lp(pos) ≈ shape(0.) rtol = 1e-14
                @test ln(-pos) ≈ shape(0.) rtol = 1e-14
            end
        end
    end

    @testset "Symmetry" begin
        for shape in test_shapes
            l = Line(0., 1., shape)
            for x in test_points
                @test l(-x) ≈ l(x) rtol = 1e-14
            end
        end
    end

    @testset "Integration" begin
        for shape in test_shapes
            l = Line(0., 1., shape)
            ref, _ = quadgk(l, -3.0, 3.0; rtol=1e-10)
            @test integrate(l, -3.0, 3.0) ≈ ref rtol = 1e-6
        end
    end
end
