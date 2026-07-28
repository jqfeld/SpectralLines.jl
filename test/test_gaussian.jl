using SpectralLines

using Test
using QuadGK

const test_sigmas = [0.01, 0.1, 1.0, 10.0, 100.0]
const test_points = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 100., 1000., 1e5]


@testset "Gaussian Line Shape" begin
    @testset "Basic functionality" begin
        g = Gaussian(1.0)

        @test g isa Gaussian{Float64}
        @test g.sigma == 1.0

        @test g(0.0) ≈ 1 / sqrt(2π) atol = 1.0e-10

        @test g(1.0) ≈ 1 / sqrt(2π) * exp(-0.5) atol = 1.0e-10
        @test g(-1.0) ≈ 1 / sqrt(2π) * exp(-0.5) atol = 1.0e-10
    end

    @testset "FWHM" begin
        for sigma in test_sigmas
            g = Gaussian(sigma)
            @test hwhm(g) ≈ fwhm(g) / 2
            @test g(fwhm(g) / 2) / g(0) ≈ 0.5 atol = 1.0e-10
        end
    end

    @testset "Normalization" begin
        for sigma in test_sigmas
            g = Gaussian(sigma)
            area, _ = quadgk(x -> g(x), -Inf, Inf, rtol=1.0e-8)
            @test area ≈ 1.0 atol = 1.0e-6
        end
    end

    @testset "Symmetry" begin
        for sigma in test_sigmas
            g = Gaussian(sigma)

            for x in test_points
                @test g(x) ≈ g(-x) atol = 1.0e-12
            end
        end
    end


    @testset "Different parameter types" begin
        g_int = Gaussian(1)
        g_float = Gaussian(1.0)
        g_rational = Gaussian(1 // 1)

        @test g_int isa Gaussian{Float64}
        @test g_float isa Gaussian{Float64}
        @test g_rational isa Gaussian{Float64}

        @test g_int(0.0) ≈ g_float(0.0) atol = 1.0e-12
        @test g_float(0.0) ≈ float(g_rational(0.0)) atol = 1.0e-12
    end

    @testset "Width scaling behavior" begin
        narrow = Gaussian(0.5)
        wide = Gaussian(2.0)

        @test narrow(0.0) > wide(0.0)

        @test narrow(0.25) > wide(0.25)
        @test narrow(1.0) < wide(1.0)
    end

    @testset "Mathematical properties" begin
        g = Gaussian(1.0)

        @test g(0.0) ≈ 1 / sqrt(2π) * exp(0) atol = 1.0e-12
        @test g(1.0) ≈ 1 / sqrt(2π) * exp(-0.5) atol = 1.0e-12
        @test g(sqrt(2)) ≈ 1 / sqrt(2π) * exp(-1.0) atol = 1.0e-12

        inflection_point = 1.0  # σ
        @test g(inflection_point) ≈ 1 / sqrt(2π) * exp(-0.5) atol = 1.0e-12
    end


    @testset "integrate" begin

        for sigma in test_sigmas
            g = Gaussian(sigma)
            # Exact erf formula — no quadrature error, so tight tolerances are valid
            @test integrate(g, -100.0 * sigma, 100.0 * sigma) ≈ 1.0 atol = 1e-10
            @test integrate(g, -2.0, 0.0) ≈ integrate(g, 0.0, 2.0) atol = 1e-12
            @test integrate(g, -5.0, 0.0) + integrate(g, 0.0, 5.0) ≈ integrate(g, -5.0, 5.0) atol = 1e-14
            for Δ in (0.01 * sigma, 0.001 * sigma)
                @test integrate(g, -Δ / 2, Δ / 2) ≈ g(0.0) * Δ rtol = 1e-4
            end
            ref, _ = quadgk(g, -3.0, 3.0; rtol=1e-10)
            @test integrate(g, -3.0, 3.0) ≈ ref rtol = 1e-6
        end
    end
end
