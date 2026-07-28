using Test, JET
using SpectralLines
import SpectralLines: faddeyeva

JET.test_package(SpectralLines; target_modules=(SpectralLines,))

const TM = (SpectralLines,)
@testset "JET call analysis" begin
    # ── Lineshapes ────────────────────────────────────────────────────────────
    @test_call target_modules = TM (Gaussian(1.0))(1.0)
    @test_call target_modules = TM fwhm(Gaussian(1.0))
    @test_call target_modules = TM integrate(Gaussian(1.0), -1.0, 1.0)



    @test_call target_modules = TM faddeyeva(1.0, 0.5)
    @test_call target_modules = TM faddeyeva(1.0 + 0.5im)
end
