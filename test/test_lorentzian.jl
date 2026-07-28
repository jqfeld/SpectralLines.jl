using SpectralLines

using Test
using QuadGK


const test_gammas = [0.01, 0.1, 1.0, 10.0, 100.0]
const test_points = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 100., 1000., 1e5]

@testset "Lorentzian Line Shape" begin
    @testset "Basic functionality" begin
        L = Lorentzian(1.0)

        @test L isa Lorentzian{Float64}
        @test L.hwhm == 1.0

        @test L(0.0) ≈ 1 / π atol = 1e-10

        @test L(1.0) ≈ 1 / π / 2 atol = 1e-10
        @test L(-1.0) ≈ 1 / π / 2 atol = 1e-10
    end

    @testset "Normalization" begin
        L = Lorentzian(1.0)
        area, _ = quadgk(x -> L(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol = 1e-6

        L2 = Lorentzian(2.0)
        area2, _ = quadgk(x -> L2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol = 1e-6
    end

    @testset "FWHM" begin
        for gamma in test_gammas
            l = Lorentzian(gamma)
            @test l(fwhm(l) / 2) / l(0) ≈ 0.5 atol = 1e-10
            @test hwhm(l) ≈ gamma
            @test hwhm(l) ≈ fwhm(l) / 2
        end
    end

    @testset "Symmetry" begin
        L = Lorentzian(1.5)

        for x in test_points
            @test L(x) ≈ L(-x) atol = 1e-12
        end
    end

    @testset "Different parameter types" begin
        L_int = Lorentzian(1)
        L_float = Lorentzian(1.0)
        L_rational = Lorentzian(1 // 1)

        @test L_int isa Lorentzian{Int}
        @test L_float isa Lorentzian{Float64}
        @test L_rational isa Lorentzian{Rational{Int}}

        @test L_int(0.0) ≈ L_float(0.0) atol = 1e-12
        @test L_float(0.0) ≈ float(L_rational(0.0)) atol = 1e-12
    end

    @testset "Width scaling behavior" begin
        narrow = Lorentzian(0.5)
        wide = Lorentzian(2.0)

        @test narrow(0.0) > wide(0.0)

        @test narrow(0.25) > wide(0.25)
        @test narrow(2.0) < wide(2.0)
    end

    @testset "Mathematical properties" begin
        L = Lorentzian(1.0)

        @test L(0.0) ≈ 1 / π atol = 1e-12
        @test L(1.0) ≈ 1 / π / 2 atol = 1e-12
        @test L(2.0) ≈ 1 / π / 5 atol = 1e-12

        hwhm = 1.0
        @test L(hwhm) ≈ 1 / π / 2 atol = 1e-12
    end

    @testset "Heavy tail behavior" begin
        L = Lorentzian(1.0)

        @test L(10.0) ≈ 1 / π / 101 atol = 1e-12
        @test L(100.0) ≈ 1 / π / 10001 atol = 1e-12

        @test L(10.0) > L(100.0)
    end

    @testset "hwhm" begin end

    @testset "integrate" begin
        for gamma in test_gammas
            l = Lorentzian(gamma)
            # Exact arctan formula — no quadrature error, so tight tolerances are valid
            @test integrate(l, -1e6 * gamma, 1e6 * gamma) ≈ 1.0 atol = 1e-4
            @test integrate(l, -2.0 * gamma, 0.0) ≈ integrate(l, 0.0, 2.0 * gamma) atol = 1e-12
            @test integrate(l, -5.0 * gamma, 0.0) + integrate(l, 0.0, 5.0 * gamma) ≈ integrate(l, -5.0 * gamma, 5.0 * gamma) atol = 1e-14
            for Δ in (0.01 * gamma, 0.001 * gamma)
                @test integrate(l, -Δ / 2, Δ / 2) ≈ l(0.0) * Δ rtol = 1e-4
            end
            ref, _ = quadgk(l, -3.0 * gamma, 3.0 * gamma; rtol=1e-10)
            @test integrate(l, -3.0 * gamma, 3.0 * gamma) ≈ ref rtol = 1e-6
        end
    end
end
