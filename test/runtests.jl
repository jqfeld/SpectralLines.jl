using SpectralLines
using Test
using Aqua
using JET

@testset "SpectralLines.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SpectralLines)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(SpectralLines; target_defined_modules = true)
    end
    # Write your tests here.
end
