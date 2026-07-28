using SpectralLines
using Test
using QuadGK
using SpecialFunctions

const test_gammas = [0.01, 0.1, 1.0, 10.0, 100.0]
const test_sigmas = [0.01, 0.1, 1.0, 10.0, 100.0]
const test_points = [0.01, 0.1, 0.5, 1.0, 2.0, 3.0, 10.0, 100., 1000., 1e5]

@testset "Voigt Line Shape" begin
    @testset "Basic functionality" begin
        V = Voigt(1.0, 0.5)

        @test V isa Voigt{Float64}
        @test V.sigma == 1.0
        @test V.gamma == 0.5

        @test V(0.0) isa Real
        @test V(0.0) > 0
    end

    @testset "Normalization" begin
        V = Voigt(1.0, 0.5)
        area, _ = quadgk(x -> V(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol = 1e-6

        V2 = Voigt(2.0, 1.0)
        area2, _ = quadgk(x -> V2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol = 1e-6
    end

    @testset "FWHM" begin
        for gamma in test_gammas
            for sigma in test_sigmas
                l = Voigt(gamma, sigma)
                @test l(fwhm(l) / 2) / l(0) ≈ 0.5 rtol = 1e-3
            end
        end
    end

    @testset "Symmetry" begin
        V = Voigt(1.5, 0.8)

        for x in test_points
            @test V(x) ≈ V(-x) atol = 1e-12
        end
    end


    @testset "Different parameter types" begin
        V_mixed = Voigt(1.0, 1)
        V_float = Voigt(1.0, 1.0)

        @test V_mixed isa Voigt{Float64}
        @test V_float isa Voigt{Float64}

        @test V_mixed(0.0) ≈ V_float(0.0) atol = 1e-12
    end

    @testset "Limiting cases" begin
        # When gamma -> 0, should approach Gaussian
        V_gaussian_like = Voigt(1.0, 1e-6)
        G = Gaussian(1.0)

        @test V_gaussian_like(0.0) ≈ G(0.0) rtol = 1e-3
        @test V_gaussian_like(1.0) ≈ G(1.0) rtol = 1e-3

        # When sigma -> 0, should approach Lorentzian
        V_lorentzian_like = Voigt(1e-6, 1.0)
        L = Lorentzian(1.0)

        @test V_lorentzian_like(0.0) ≈ L(0.0) rtol = 1e-3
        @test V_lorentzian_like(1.0) ≈ L(1.0) rtol = 1e-3
    end

    @testset "Mathematical properties" begin
        V = Voigt(1.0, 0.5)

        @test V(0.0) > V(1.0)
        @test V(1.0) > V(2.0)

        @test all(V(x) > 0 for x in [-5, -1, 0, 1, 5])

        center_val = V(0.0)
        @test center_val isa Real
        @test center_val > 0
    end

    @testset "Faddeeva function usage" begin
        # Test that the Voigt profile uses the Faddeeva function correctly

        for σ in test_sigmas, γ in test_gammas
            for x in test_points
                V = Voigt(σ, γ)
                expected = real(faddeeva((x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)

                @test V(x) ≈ expected rtol = 1e-4
            end
        end
    end

    @testset "hwhm" begin
        @test hwhm(Voigt(1.0, 0.5)) ≈ fwhm(Voigt(1.0, 0.5)) / 2
    end

    @testset "integrate" begin
        # GH branch (y = γ/(σ√2) ≥ 0.7): exact additivity and accurate for wide intervals
        v_gh = Voigt(1.0, 2.0)   # y ≈ 1.41
        @test integrate(v_gh, -5.0, 0.0) + integrate(v_gh, 0.0, 5.0) ≈ integrate(v_gh, -5.0, 5.0) atol = 1e-15
        ref_wide, _ = quadgk(v_gh, -100.0, 100.0; rtol=1e-10)
        @test integrate(v_gh, -100.0, 100.0) ≈ ref_wide rtol = 1e-4

        # GL branch (y < 0.7): Gaussian-dominated, accurate for bins ≤ 2×FWHM
        v_g = Voigt(1.0, 0.01)   # nearly Gaussian
        hw_g = hwhm(v_g)
        @test sum(integrate(v_g, (i - 1) * hw_g - hw_g / 2, (i - 1) * hw_g + hw_g / 2) for i in -50:50) ≈ 1.0 atol = 2e-3

        # Symmetry
        v = Voigt(1.0, 0.5)
        @test integrate(v, -2.0, 0.0) ≈ integrate(v, 0.0, 2.0) atol = 1e-8

        # Accuracy vs QuadGK across both branches: < 0.5% for bins up to 1×FWHM
        for (σ, γ) in [(1.0, 0.3), (1.0, 1.0), (1.0, 5.0)]
            v = Voigt(σ, γ)
            hw = hwhm(v)
            for frac in [0.2, 0.5, 1.0], ctr in [0.0, 2hw]
                D = frac * 2hw
                ref, _ = quadgk(v, ctr - D / 2, ctr + D / 2; rtol=1e-10)
                @test integrate(v, ctr - D / 2, ctr + D / 2) ≈ ref rtol = 5e-3
            end
        end

        # Lorentzian limit (σ → 0): Voigt collapses to Lorentzian
        l = Lorentzian(1.0)
        v_l = Voigt(1e-8, 1.0)
        for (a, b) in [(-2.0, 2.0), (0.5, 3.0)]
            @test integrate(v_l, a, b) ≈ integrate(l, a, b) rtol = 1e-3
        end
    end

    @testset "Exact degenerate limits (σ=0, γ=0) — the raw callable" begin
        # Before the σ==0/γ==0 special-casing in `(::Voigt)(x)`, σ=0 produced NaN
        # (inv(0) = Inf, and x=0 gives 0*Inf = NaN) and γ=0 was merely ~1e-6 off
        # from the exact Gaussian (faddeyeva's own approximation error at y=0) —
        # see _research/TODO.md. These are exact-equality checks (not `integrate`,
        # which already special-cases these via its own quadrature branches and
        # was passing before this fix) against the raw `Voigt(...)(x)` callable.
        xs = [-5.0, -2.0, -0.5, 0.0, 0.3, 1.0, 4.0]

        @testset "σ = 0 → exact Lorentzian" begin
            γ = 0.35
            v = Voigt(0.0, γ)
            l = Lorentzian(γ)
            @test all(isfinite, v.(xs))
            for x in xs
                @test v(x) === l(x)   # bit-identical: both dispatch to the same `lorentzian(x, γ)`
            end
        end

        @testset "γ = 0 → exact Gaussian" begin
            σ = 0.4
            v = Voigt(σ, 0.0)
            g = Gaussian(σ)
            @test all(isfinite, v.(xs))
            for x in xs
                @test v(x) === g(x)   # bit-identical: both dispatch to the same `gaussian(x, σ)`
            end
        end

        @testset "σ = 0 and γ = 0 (doubly-degenerate, σ branch checked first)" begin
            # Voigt(0,0) takes the sigma==0 branch and reduces to lorentzian(x, 0.0),
            # inheriting whatever Lorentzian(0.0) itself already does at these points
            # (pre-existing behaviour, not something this fix changes or should try
            # to improve): NaN at the centre (0/0), exactly 0 elsewhere (0/x²).
            v = Voigt(0.0, 0.0)
            @test isnan(v(0.0)) && isnan(Lorentzian(0.0)(0.0))
            @test v(1.0) === Lorentzian(0.0)(1.0) === 0.0
        end
    end
end
