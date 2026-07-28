using Test
using SpectralLines
import SpectralLines: faddeyeva

# For each region, verify faddeyeva(x, y) == real(faddeyeva(x + im*y)).
# This checks the real-only helpers against the complex reference implementation.
# Region boundaries (s = x²+y²):
#   1: s >= 38000   2: s >= 256   3: s >= 62
#   4: s >= 30, y² >= 1e-13
#   5a: s > 2.5, y² < 5e-9
#   5b: s > 2.5, 5e-9 <= y² < 0.072
#   6: everything else

@testset "faddeyeva real-only vs complex reference" begin

    tol = 1e-12

    @testset "region 1 (s >= 38000)" begin
        for (x, y) in [(200.0, 0.1), (0.0, 196.0), (150.0, 100.0)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 2 (256 <= s < 38000)" begin
        for (x, y) in [(16.0, 0.5), (10.0, 10.0), (0.0, 16.1)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 3 (62 <= s < 256)" begin
        for (x, y) in [(8.0, 0.1), (5.0, 6.0), (0.0, 8.0)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 4 (30 <= s < 62, y² >= 1e-13)" begin
        for (x, y) in [(5.5, 0.5), (5.0, 1.0), (0.0, 5.48)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 5a (s > 2.5, y² < 5e-9)" begin
        for (x, y) in [(2.0, 1e-5), (3.0, 1e-6), (1.7, 2e-5)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 5b (s > 2.5, 5e-9 <= y² < 0.072)" begin
        for (x, y) in [(2.0, 0.05), (1.8, 0.1), (3.0, 0.01)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

    @testset "region 6 (s <= 2.5 or y² >= 0.072)" begin
        for (x, y) in [(0.5, 0.3), (1.0, 0.5), (0.0, 1.0), (1.2, 0.8)]
            @test faddeyeva(x, y) ≈ real(faddeyeva(x + im*y)) atol=tol
        end
    end

end
