using SpectralLines
using Test
using QuadGK

# const test_gammas = [0.01, 0.1, 1.0, 10.0, 100.0]
# const test_sigmas = [0.01, 0.1, 1.0, 10.0, 100.0]
const test_points = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 100., 1000., 1e5]

@testset "PseudoVoigt line shape" begin
    @testset "Basic functionality" begin
        VA = PseudoVoigt(1.0, 0.5)

        @test VA isa PseudoVoigt
        @test VA.sigma == 1.0
        @test VA.gamma == 0.5

        @test VA(0.0) isa Real
        @test VA(0.0) > 0
    end

    @testset "Normalization" begin
        VA = PseudoVoigt(1.0, 0.5)
        area, _ = quadgk(x -> VA(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol = 1e-5

        VA2 = PseudoVoigt(2.0, 1.0)
        area2, _ = quadgk(x -> VA2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol = 1e-5
    end

    @testset "Symmetry" begin
        VA = PseudoVoigt(1.5, 0.8)

        for x in test_points
            @test VA(x) ≈ VA(-x) atol = 1e-12
        end
    end

    @testset "Different parameter types" begin
        VA_mixed = PseudoVoigt(1.0, 1)
        VA_float = PseudoVoigt(1.0, 1.0)

        @test VA_mixed(0.0) ≈ VA_float(0.0) atol = 1e-12
    end

    @testset "Limiting cases" begin
        # When gamma -> 0, should approach Gaussian-like behavior
        VA_gaussian_like = PseudoVoigt(1.0, 1e-6)
        G = Gaussian(1.0)

        @test VA_gaussian_like(0.0) ≈ G(0.0) rtol = 1e-2

        # When sigma -> 0, should approach Lorentzian-like behavior
        VA_lorentzian_like = PseudoVoigt(1e-6, 1.0)
        L = Lorentzian(1.0)

        @test VA_lorentzian_like(0.0) ≈ L(0.0) rtol = 1e-2
    end

    @testset "Approximation accuracy" begin
        # Test that PseudoVoigt is reasonably close to exact Voigt
        σ = 1.0
        γ = 0.5

        VA = PseudoVoigt(σ, γ)
        V = Voigt(σ, γ)

        test_points = [-2.0, -1.0, 0.0, 1.0, 2.0]

        for x in test_points
            @test VA(x) ≈ V(x) rtol = 2e-2  # Allow 2% relative error
        end
    end

    @testset "Mathematical properties" begin
        VA = PseudoVoigt(1.0, 0.5)

        @test VA(0.0) > VA(1.0)
        @test VA(1.0) > VA(2.0)

        @test all(VA(x) > 0 for x in [-5, -1, 0, 1, 5])

        center_val = VA(0.0)
        @test center_val isa Real
        @test center_val > 0
    end

    @testset "Width parameter conversion" begin
        # Test that the approximation handles FWHM conversion correctly
        σ = 1.0
        γ = 0.5

        VA = PseudoVoigt(σ, γ)

        # The approximation should produce reasonable values
        @test VA(0.0) > 0
        @test VA(1.0) > 0
        @test VA(2.0) > 0

        # Check monotonic decrease from center
        @test VA(0.0) > VA(0.5) > VA(1.0)
    end

    @testset "Edge cases" begin
        # Test with very small parameters
        VA_small = PseudoVoigt(1e-3, 1e-3)
        @test VA_small(0.0) > 0

        # Test with large parameters
        VA_large = PseudoVoigt(10.0, 5.0)
        @test VA_large(0.0) > 0

        # Test with unequal parameters
        VA_unequal = PseudoVoigt(0.1, 5.0)
        @test VA_unequal(0.0) > 0
    end

    @testset "Mixing models" begin
        σ, γ = 1.0, 0.5
        V = Voigt(σ, γ)
        models = [Thompson(), Kielkopf(), Olivero(), LiuLin(), Rocco()]

        for m in models
            VA = PseudoVoigt(σ, γ, m)
            # Each model must be callable and return a positive real
            @test VA(0.0) isa Real
            @test VA(0.0) > 0
            # Must be symmetric
            @test VA(1.0) ≈ VA(-1.0) atol = 1e-12
            # Must be normalised within 1%
            area, _ = quadgk(x -> VA(x), -Inf, Inf, rtol=1e-8)
            @test area ≈ 1.0 atol = 1e-2
            # Must approximate the exact Voigt within 5%
            for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
                @test VA(x) ≈ V(x) rtol = 5e-2
            end
            # fwhm must be positive
            @test fwhm(VA) > 0
        end

        # Default constructor still gives Thompson
        @test PseudoVoigt(σ, γ)(0.0) ≈ PseudoVoigt(σ, γ, Thompson())(0.0) atol = 1e-15

        # Rocco large-y limit: pure Lorentzian (η → 1)
        VA_rocco_large_y = PseudoVoigt(1e-4, 1.0, Rocco())
        L = Lorentzian(1.0)
        @test VA_rocco_large_y(0.0) ≈ L(0.0) rtol = 1e-2
    end


    @testset "integrate" begin
        # Analytical formula (erf + arctan combination) — exact for any bin width
        for M in (Thompson(), Kielkopf(), Olivero(), LiuLin(), Rocco())
            pv = PseudoVoigt(1.0, 0.5, M)
            @test integrate(pv, -1000.0, 1000.0) ≈ 1.0 atol = 1e-3
            @test integrate(pv, -2.0, 0.0) ≈ integrate(pv, 0.0, 2.0) atol = 1e-8
            @test integrate(pv, -5.0, 0.0) + integrate(pv, 0.0, 5.0) ≈ integrate(pv, -5.0, 5.0) atol = 1e-14
            for Δ in (0.01, 0.001)
                @test integrate(pv, -Δ / 2, Δ / 2) ≈ pv(0.0) * Δ rtol = 1e-3
            end
        end
    end
end
